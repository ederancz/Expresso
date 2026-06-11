"""Phase 4 — outlier QC workbook for intrinsic master outputs."""

from __future__ import annotations

import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from openpyxl import Workbook

from .area import normalize_region
from .metadata import is_exclude_flag
from .values import coerce_param_value

QC_WORKBOOK_NAME = "Intrinsic_QC.xlsx"

STRATIFIED_GROUPS = ("VISp-L5", "V2M-L5", "VISp-L2/3", "V2M-L2/3")

_IQR_MULTIPLIER = 3.0

_LOG_SKEW_RE = re.compile(r"(?:^|_)(rin|rheobase|tau)(?:_|$)", re.IGNORECASE)

GROUP_MEAN_ID = "__GROUP_MEAN__"
GROUP_SD_ID = "__GROUP_SD__"

# Documented in README.md and METHODS.md — keep in sync.
FLAG_TYPE_DEFINITIONS: dict[str, str] = {
    "iqr": (
        "Value outside the within-group 3×IQR fence (Q1 − 3×IQR or Q3 + 3×IQR). "
        "Group = region–layer (control) or region–layer|experiment (pharmacology)."
    ),
    "iqr+global": (
        "Within-group 3×IQR outlier that is also outside the pooled global 3×IQR fence "
        "for that parameter across all groups in the same scope."
    ),
    "log_iqr": (
        "On log10-transformed values (Rin, rheobase, τ only; positive values), "
        "outside the within-group 3×IQR fence on the log scale."
    ),
    "log_iqr+global": "Log-scale 3×IQR outlier that is also outside the pooled global log-scale fence.",
    "bio:Rin<0": "Input resistance below zero.",
    "bio:tau<0": "Membrane time constant below zero.",
    "bio:AP_peak<-20mV": "Action-potential peak more depolarized than −20 mV.",
    "bio:Vm_rest>-30mV": "Resting membrane potential above −30 mV.",
    "bio:rheobase_out_of_range": "Rheobase below 0 pA or above 2000 pA.",
    "bio:AP_width<0.1ms": "AP width below 0.1 ms.",
    "bio:sag_pct_outside_(0,100)": "Sag percentage not strictly between 0 and 100.",
    "bio:ratio_outside_(0,1)": "Adaptation or burst ratio not strictly between 0 and 1.",
    "bio:res_freq_out_of_range": "Resonance frequency below 0 Hz or above 200 Hz.",
}


def _value_col(param: str) -> str:
    return f"{param}__value"


def _flag_matrix_columns(param_columns: list[str]) -> list[str]:
    cols: list[str] = []
    for param in param_columns:
        cols.append(param)
        cols.append(_value_col(param))
    return cols


@dataclass
class DistributionStats:
    n: int = 0
    mean: float | None = None
    sd: float | None = None
    median: float | None = None
    q1: float | None = None
    q3: float | None = None
    min: float | None = None
    max: float | None = None
    iqr: float | None = None
    lower_fence: float | None = None
    upper_fence: float | None = None
    log_lower_fence: float | None = None
    log_upper_fence: float | None = None
    log_n: int = 0


@dataclass
class CellFlags:
    cell_id: str
    region: str = ""
    layer: str = ""
    group: str = ""
    scope: str = "control"
    experiment: str = ""
    total_flags: int = 0
    param_flags: dict[str, str] = field(default_factory=dict)
    param_values: dict[str, float] = field(default_factory=dict)


def _group_key(region: Any, layer: Any) -> str:
    reg = normalize_region(str(region or "").strip())
    lay = str(layer or "").strip()
    if lay == "L2-3":
        lay = "L2/3"
    if reg in {"VISp", "V2M"} and lay in {"L5", "L2/3"}:
        return f"{reg}-{lay}"
    return ""


def _param_tail(param: str) -> str:
    if "__" in param:
        return param.split("__", 1)[1]
    return param


def _param_tail_lower(param: str) -> str:
    return _param_tail(param).lower()


def is_log_skewed_param(param: str) -> bool:
    return bool(_LOG_SKEW_RE.search(_param_tail_lower(param)))


def _distribution_stats(values: list[float], *, log_scale: bool = False) -> DistributionStats:
    stats = DistributionStats()
    if not values:
        return stats

    arr = np.asarray(values, dtype=float)
    stats.n = int(arr.size)
    stats.mean = float(np.mean(arr))
    stats.sd = float(np.std(arr, ddof=1)) if stats.n > 1 else 0.0
    stats.median = float(np.median(arr))
    stats.q1 = float(np.percentile(arr, 25))
    stats.q3 = float(np.percentile(arr, 75))
    stats.min = float(np.min(arr))
    stats.max = float(np.max(arr))
    stats.iqr = stats.q3 - stats.q1
    if stats.iqr == 0:
        stats.lower_fence = stats.q1
        stats.upper_fence = stats.q3
    else:
        stats.lower_fence = stats.q1 - _IQR_MULTIPLIER * stats.iqr
        stats.upper_fence = stats.q3 + _IQR_MULTIPLIER * stats.iqr

    if log_scale:
        positive = arr[arr > 0]
        stats.log_n = int(positive.size)
        if stats.log_n >= 2:
            log_arr = np.log10(positive)
            lq1 = float(np.percentile(log_arr, 25))
            lq3 = float(np.percentile(log_arr, 75))
            liqr = lq3 - lq1
            if liqr == 0:
                stats.log_lower_fence = lq1
                stats.log_upper_fence = lq3
            else:
                stats.log_lower_fence = lq1 - _IQR_MULTIPLIER * liqr
                stats.log_upper_fence = lq3 + _IQR_MULTIPLIER * liqr
    return stats


def biological_limit_violation(param: str, value: float) -> str | None:
    """Return a short reason string when value violates plausibility limits."""
    tail = _param_tail_lower(param)

    if "rin" in tail and value < 0:
        return "Rin<0"

    if re.search(r"(?:^|_)tau(?:_|$)", tail) and value < 0:
        return "tau<0"

    if "ap_peak" in tail and value < -20:
        return "AP_peak<-20mV"

    if "vm_rest" in tail and value > -30:
        return "Vm_rest>-30mV"

    if "rheobase" in tail and (value < 0 or value > 2000):
        return "rheobase_out_of_range"

    if "ap_width" in tail and value < 0.1:
        return "AP_width<0.1ms"

    if "sag_percentage" in tail and not (0 < value < 100):
        return "sag_pct_outside_(0,100)"

    if "adaptation_ratio" in tail or "burst_ratio" in tail:
        if not (0 < value < 1):
            return "ratio_outside_(0,1)"

    if ("res._freq" in tail or "resonance" in tail) and (value < 0 or value > 200):
        return "res_freq_out_of_range"

    return None


def _is_iqr_outlier(value: float, stats: DistributionStats) -> bool:
    if stats.n < 4 or stats.lower_fence is None or stats.upper_fence is None:
        return False
    return value < stats.lower_fence or value > stats.upper_fence


def _is_log_iqr_outlier(value: float, stats: DistributionStats) -> bool:
    if (
        stats.log_n < 4
        or stats.log_lower_fence is None
        or stats.log_upper_fence is None
        or value <= 0
    ):
        return False
    log_val = math.log10(value)
    return log_val < stats.log_lower_fence or log_val > stats.log_upper_fence


def _format_flag_parts(parts: list[str]) -> str:
    return ";".join(parts) if parts else ""


def _collect_numeric_by_group(
    rows: list[dict[str, Any]],
    param_columns: list[str],
    *,
    group_fn,
    exclude_from_distribution: bool = True,
) -> tuple[dict[str, dict[str, list[tuple[str, float]]]], dict[str, dict[str, list[float]]]]:
    """
    Return (all_values, distribution_values).

    all_values: group -> param -> [(row_key, value), ...] for every parseable cell.
    distribution_values: group -> param -> [value, ...] excluding flagged rows when requested.
    """
    all_values: dict[str, dict[str, list[tuple[str, float]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    dist_values: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))

    for row in rows:
        group = group_fn(row)
        if not group:
            continue
        row_key = row.get("cell_id", "")
        excluded = exclude_from_distribution and is_exclude_flag(row.get("exclude_flag"))
        for param in param_columns:
            val = coerce_param_value(row.get(param))
            if val is None:
                continue
            all_values[group][param].append((row_key, val))
            if not excluded:
                dist_values[group][param].append(val)
    return all_values, dist_values


def _compute_group_param_stats(
    dist_values: dict[str, dict[str, list[float]]],
    param_columns: list[str],
    *,
    groups: tuple[str, ...],
) -> dict[tuple[str, str], DistributionStats]:
    out: dict[tuple[str, str], DistributionStats] = {}
    for group in groups:
        for param in param_columns:
            vals = dist_values.get(group, {}).get(param, [])
            out[(group, param)] = _distribution_stats(vals, log_scale=is_log_skewed_param(param))
    return out


def _compute_global_param_stats(
    dist_values: dict[str, dict[str, list[float]]],
    param_columns: list[str],
) -> dict[str, DistributionStats]:
    pooled: dict[str, list[float]] = defaultdict(list)
    for group_vals in dist_values.values():
        for param, vals in group_vals.items():
            pooled[param].extend(vals)
    return {
        param: _distribution_stats(pooled.get(param, []), log_scale=is_log_skewed_param(param))
        for param in param_columns
    }


def _pharma_group_key(row: dict[str, Any]) -> str:
    base = _group_key(row.get("region"), row.get("layer"))
    experiment = str(row.get("experiment", "") or "").strip()
    if not base:
        return ""
    if experiment:
        return f"{base}|{experiment}"
    return base


def _build_param_summary_rows(
    *,
    scope: str,
    groups: list[str],
    param_columns: list[str],
    group_stats: dict[tuple[str, str], DistributionStats],
    flagged_cells: list[CellFlags],
) -> list[dict[str, Any]]:
    counts: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for cell in flagged_cells:
        if cell.scope != scope:
            continue
        for param, flag_str in cell.param_flags.items():
            key = (cell.group, param)
            for part in flag_str.split(";"):
                if part == "iqr":
                    counts[key]["n_iqr"] += 1
                elif part == "iqr+global":
                    counts[key]["n_iqr"] += 1
                    counts[key]["n_also_global"] += 1
                elif part == "log_iqr":
                    counts[key]["n_log_iqr"] += 1
                elif part == "log_iqr+global":
                    counts[key]["n_log_iqr"] += 1
                    counts[key]["n_also_global"] += 1
                elif part.startswith("bio:"):
                    counts[key]["n_bio"] += 1

    rows: list[dict[str, Any]] = []
    for group in groups:
        for param in param_columns:
            stats = group_stats.get((group, param), DistributionStats())
            c = counts.get((group, param), {})
            rows.append(
                {
                    "scope": scope,
                    "group": group,
                    "parameter": param,
                    "n": stats.n,
                    "mean": stats.mean,
                    "SD": stats.sd,
                    "median": stats.median,
                    "Q1": stats.q1,
                    "Q3": stats.q3,
                    "min": stats.min,
                    "max": stats.max,
                    "IQR": stats.iqr,
                    "lower_fence": stats.lower_fence,
                    "upper_fence": stats.upper_fence,
                    "log_n": stats.log_n,
                    "log_lower_fence": stats.log_lower_fence,
                    "log_upper_fence": stats.log_upper_fence,
                    "n_iqr_outliers": c.get("n_iqr", 0),
                    "n_log_iqr_outliers": c.get("n_log_iqr", 0),
                    "n_bio_outliers": c.get("n_bio", 0),
                    "n_also_global_outliers": c.get("n_also_global", 0),
                }
            )
    return rows


def _rows_to_cell_flags(
    rows: list[dict[str, Any]],
    param_columns: list[str],
    *,
    scope: str,
    group_fn,
    group_stats: dict[tuple[str, str], DistributionStats],
    global_stats: dict[str, DistributionStats],
) -> list[CellFlags]:
    """Build a flag entry for every row in scope (empty flags when none)."""
    out: list[CellFlags] = []
    for row in rows:
        group = group_fn(row)
        if not group:
            continue
        cell = CellFlags(
            cell_id=str(row.get("cell_id", "")),
            region=normalize_region(str(row.get("region", ""))),
            layer=str(row.get("layer", "")),
            group=group,
            scope=scope,
            experiment=str(row.get("experiment", "") or ""),
        )
        for param in param_columns:
            val = coerce_param_value(row.get(param))
            if val is None:
                continue
            parts: list[str] = []
            gstats = group_stats.get((group, param), DistributionStats())
            if _is_iqr_outlier(val, gstats):
                if _is_iqr_outlier(val, global_stats.get(param, DistributionStats())):
                    parts.append("iqr+global")
                else:
                    parts.append("iqr")
            if is_log_skewed_param(param):
                if _is_log_iqr_outlier(val, gstats):
                    if _is_log_iqr_outlier(val, global_stats.get(param, DistributionStats())):
                        parts.append("log_iqr+global")
                    else:
                        parts.append("log_iqr")
            bio = biological_limit_violation(param, val)
            if bio:
                parts.append(f"bio:{bio}")
            if parts:
                cell.param_flags[param] = _format_flag_parts(parts)
                cell.param_values[param] = val
                cell.total_flags += len(parts)
        out.append(cell)
    return out


def _flagged_only(cells: list[CellFlags]) -> list[CellFlags]:
    return [c for c in cells if c.param_flags]


def _parse_group_parts(group: str) -> tuple[str, str, str]:
    """Return (region, layer, experiment) from a group key."""
    if "|" in group:
        base, experiment = group.split("|", 1)
    else:
        base, experiment = group, ""
    if "-" in base:
        region, layer = base.split("-", 1)
    else:
        region, layer = base, ""
    return region, layer, experiment


def _build_group_summary_rows(
    *,
    scope: str,
    groups: list[str],
    param_columns: list[str],
    group_stats: dict[tuple[str, str], DistributionStats],
) -> list[dict[str, Any]]:
    """Two header rows per group: group mean and group SD in the __value twin columns."""
    rows: list[dict[str, Any]] = []
    for group in groups:
        region, layer, experiment = _parse_group_parts(group)
        for row_id, stat_attr in ((GROUP_MEAN_ID, "mean"), (GROUP_SD_ID, "sd")):
            row: dict[str, Any] = {
                "cell_id": row_id,
                "scope": scope,
                "region": region,
                "layer": layer,
                "group": group,
                "experiment": experiment,
                "total_flags": "",
            }
            for param in param_columns:
                row[param] = ""
                stats = group_stats.get((group, param), DistributionStats())
                row[_value_col(param)] = getattr(stats, stat_attr, None) if stats.n else ""
            rows.append(row)
    return rows


def _build_flag_matrix_rows(cells: list[CellFlags], param_columns: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for cell in cells:
        row: dict[str, Any] = {
            "cell_id": cell.cell_id,
            "scope": cell.scope,
            "region": cell.region,
            "layer": cell.layer,
            "group": cell.group,
            "experiment": cell.experiment,
            "total_flags": cell.total_flags,
        }
        for param in param_columns:
            row[param] = cell.param_flags.get(param, "")
            row[_value_col(param)] = cell.param_values.get(param, "") if param in cell.param_flags else ""
        rows.append(row)
    return rows


def _build_suspicious_rows(flagged: list[CellFlags]) -> list[dict[str, Any]]:
    ranked = sorted(flagged, key=lambda c: (-c.total_flags, c.cell_id, c.scope))
    rows: list[dict[str, Any]] = []
    for rank, cell in enumerate(ranked, start=1):
        flagged_params = [
            f"{param}:{flags}" for param, flags in sorted(cell.param_flags.items())
        ]
        rows.append(
            {
                "rank": rank,
                "cell_id": cell.cell_id,
                "scope": cell.scope,
                "region": cell.region,
                "layer": cell.layer,
                "group": cell.group,
                "experiment": cell.experiment,
                "total_flags": cell.total_flags,
                "flagged_parameters": " | ".join(flagged_params),
            }
        )
    return rows


def _write_qc_workbook(
    path: Path,
    *,
    flag_matrix_rows: list[dict[str, Any]],
    param_summary_rows: list[dict[str, Any]],
    suspicious_rows: list[dict[str, Any]],
    flag_matrix_columns: list[str],
    flag_matrix_freeze_row: int,
) -> None:
    wb = Workbook()
    wb.remove(wb.active)

    summary_cols = [
        "scope",
        "group",
        "parameter",
        "n",
        "mean",
        "SD",
        "median",
        "Q1",
        "Q3",
        "min",
        "max",
        "IQR",
        "lower_fence",
        "upper_fence",
        "log_n",
        "log_lower_fence",
        "log_upper_fence",
        "n_iqr_outliers",
        "n_log_iqr_outliers",
        "n_bio_outliers",
        "n_also_global_outliers",
    ]
    suspicious_cols = [
        "rank",
        "cell_id",
        "scope",
        "region",
        "layer",
        "group",
        "experiment",
        "total_flags",
        "flagged_parameters",
    ]

    for title, rows, columns, freeze_row in (
        ("flag_matrix", flag_matrix_rows, flag_matrix_columns, flag_matrix_freeze_row),
        ("param_summary", param_summary_rows, summary_cols, 2),
        ("suspicious_cells", suspicious_rows, suspicious_cols, 2),
    ):
        ws = wb.create_sheet(title)
        ws.append(columns)
        ws.freeze_panes = f"A{freeze_row}"
        for row in rows:
            ws.append([row.get(c) for c in columns])

    wb.save(path)


def write_qc_workbook(
    output_dir: Path,
    control_rows: list[dict[str, Any]],
    effect_rows: list[dict[str, Any]],
    meta_columns: list[str],
    param_columns: list[str],
) -> dict[str, Any]:
    """
    Run Phase 4 outlier QC and write ``Intrinsic_QC.xlsx``.

    Returns a summary dict suitable for ``run_manifest``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    qc_path = output_dir / QC_WORKBOOK_NAME

    stratified = STRATIFIED_GROUPS

    _, control_dist = _collect_numeric_by_group(
        control_rows,
        param_columns,
        group_fn=lambda r: _group_key(r.get("region"), r.get("layer")),
    )
    control_group_stats = _compute_group_param_stats(
        control_dist, param_columns, groups=stratified
    )
    control_global_stats = _compute_global_param_stats(control_dist, param_columns)

    control_group_fn = lambda r: _group_key(r.get("region"), r.get("layer"))

    control_cells = _rows_to_cell_flags(
        control_rows,
        param_columns,
        scope="control",
        group_fn=control_group_fn,
        group_stats=control_group_stats,
        global_stats=control_global_stats,
    )
    control_flagged = _flagged_only(control_cells)

    pharma_groups = sorted(
        {
            _pharma_group_key(r)
            for r in effect_rows
            if _pharma_group_key(r)
        }
    )
    _, effect_dist = _collect_numeric_by_group(
        effect_rows,
        param_columns,
        group_fn=_pharma_group_key,
    )
    effect_group_stats = _compute_group_param_stats(
        effect_dist, param_columns, groups=tuple(pharma_groups)
    )
    effect_global_stats = _compute_global_param_stats(effect_dist, param_columns)

    effect_cells = _rows_to_cell_flags(
        effect_rows,
        param_columns,
        scope="pharmacology",
        group_fn=_pharma_group_key,
        group_stats=effect_group_stats,
        global_stats=effect_global_stats,
    )
    effect_flagged = _flagged_only(effect_cells)

    param_summary_rows = _build_param_summary_rows(
        scope="control",
        groups=list(stratified),
        param_columns=param_columns,
        group_stats=control_group_stats,
        flagged_cells=control_flagged,
    )
    param_summary_rows.extend(
        _build_param_summary_rows(
            scope="pharmacology",
            groups=pharma_groups,
            param_columns=param_columns,
            group_stats=effect_group_stats,
            flagged_cells=effect_flagged,
        )
    )

    all_flagged = control_flagged + effect_flagged
    flag_matrix_meta = [
        "cell_id",
        "scope",
        "region",
        "layer",
        "group",
        "experiment",
        "total_flags",
    ]
    flag_param_cols = _flag_matrix_columns(param_columns)
    flag_matrix_columns = flag_matrix_meta + flag_param_cols

    summary_rows: list[dict[str, Any]] = []
    summary_rows.extend(
        _build_group_summary_rows(
            scope="control",
            groups=list(stratified),
            param_columns=param_columns,
            group_stats=control_group_stats,
        )
    )
    summary_rows.extend(
        _build_group_summary_rows(
            scope="pharmacology",
            groups=pharma_groups,
            param_columns=param_columns,
            group_stats=effect_group_stats,
        )
    )
    flag_matrix_rows = summary_rows + _build_flag_matrix_rows(
        control_cells + effect_cells, param_columns
    )
    suspicious_rows = _build_suspicious_rows(all_flagged)
    flag_matrix_freeze_row = 2 + len(summary_rows)

    _write_qc_workbook(
        qc_path,
        flag_matrix_rows=flag_matrix_rows,
        param_summary_rows=param_summary_rows,
        suspicious_rows=suspicious_rows,
        flag_matrix_columns=flag_matrix_columns,
        flag_matrix_freeze_row=flag_matrix_freeze_row,
    )

    n_excluded_from_dist = sum(
        1
        for r in control_rows + effect_rows
        if is_exclude_flag(r.get("exclude_flag"))
    )
    n_control_flagged = len(control_flagged)
    n_effect_flagged = len(effect_flagged)
    n_total_flags = sum(c.total_flags for c in all_flagged)

    report: dict[str, Any] = {
        "qc_workbook": str(qc_path.resolve()),
        "stratified_groups": list(stratified),
        "pharmacology_groups": pharma_groups,
        "parameters_tested": len(param_columns),
        "cells_excluded_from_iqr_distributions": n_excluded_from_dist,
        "flagged_cells": {
            "control": n_control_flagged,
            "pharmacology": n_effect_flagged,
            "total": n_control_flagged + n_effect_flagged,
        },
        "total_flag_events": n_total_flags,
        "flag_matrix_summary_rows": len(summary_rows),
        "flag_matrix_value_twin_columns": True,
        "flag_types": FLAG_TYPE_DEFINITIONS,
        "sheets": ["flag_matrix", "param_summary", "suspicious_cells"],
        "methods": {
            "iqr_fence": "Q1 - 3×IQR / Q3 + 3×IQR (within-group; exclude_flag=1 omitted from distributions)",
            "log_iqr_params": "Rin, rheobase, tau (log10 scale, positive values only)",
            "global_note": "within-group IQR outliers also outside global pooled IQR are tagged +global",
            "pharmacology": "absolute effect values; groups = region-layer|experiment",
            "flag_matrix_layout": (
                "Top rows: __GROUP_MEAN__ / __GROUP_SD__ per scope×group with twin __value columns; "
                "then cell rows with flag codes and flagged values in adjacent __value columns."
            ),
        },
    }
    return report
