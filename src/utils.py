"""Gene resolution, brain-area mapping, and helper utilities."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

# Cortical dissection ROI acronyms (WMB-10X region_of_interest_acronym).
_CORTICAL_ROIS = frozenset({
    "ACA", "AI", "AUD", "AUD-TEa-PERI-ECT", "MO-FRP", "MOp", "PL-ILA-ORB",
    "RSP", "SS-GU-VISC", "SSp", "TEa-PERI-ECT", "VIS", "VIS-PTLp",
})

# CCF v3 sub-regions -> WMB-10X dissection ROI (scRNA lacks per-cell CCF labels).
_CCF_TO_DISSECTION_ROI: dict[str, str] = {
    "VISp": "VIS",
    "VISpm": "VIS",
    "VISam": "VIS",
    "VISl": "VIS",
    "VISal": "VIS",
    "VISrl": "VIS",
    "VISli": "VIS",
    "RSPagl": "RSP",
    "MOs": "MO-FRP",
    "CP": "STRd",
    "ACB": "STRv",
    "CA1": "HIP",
    "DG": "HIP",
}

ALLEN_CCF_DIR = "Allen-CCF-2020"

# Parcellation columns to try (finest first) when assigning MERFISH brain_area.
_MERFISH_PARCELLATION_LEVELS = (
    "parcellation_substructure",
    "parcellation_structure",
    "parcellation_division",
    "parcellation_category",
)

# feature_matrix_label suffix -> config brain_area (WMB-10Xv3 packages).
_PACKAGE_TO_BRAIN_AREA = {
    "CB": "CB",
    "CTXsp": "CTX",
    "HPF": "HIP",
    "HY": "HY",
    "Isocortex-1": "CTX",
    "Isocortex-2": "CTX",
    "MB": "MB",
    "MY": "MY",
    "OLF": "OLF",
    "P": "P",
    "PAL": "PAL",
    "STR": "STR",
    "TH": "TH",
}


def assign_merfish_brain_area(
    cell_meta: pd.DataFrame,
    brain_areas: list[str],
) -> pd.Series:
    """
    Map each MERFISH cell to a config brain_area using CCF parcellation columns.

    Uses the finest parcellation level that matches any requested area.
    """
    area_set = set(brain_areas)
    assigned = pd.Series(np.nan, index=cell_meta.index, dtype=object)
    for col in _MERFISH_PARCELLATION_LEVELS:
        if col not in cell_meta.columns:
            continue
        values = cell_meta[col].astype(str)
        mask = values.isin(area_set) & assigned.isna()
        assigned[mask] = values[mask]

    if assigned.notna().any():
        return assigned

    warnings.warn(
        "Could not assign MERFISH brain_area from parcellation columns; "
        f"expected one of {_MERFISH_PARCELLATION_LEVELS}",
        UserWarning,
        stacklevel=2,
    )
    return pd.Series(np.nan, index=cell_meta.index, dtype=object)


def resolve_gene_ids(
    gene_df: pd.DataFrame,
    symbols: list[str],
    symbol_col: str = "gene_symbol",
) -> tuple[dict[str, str], list[str]]:
    """
    Map gene symbols to Ensembl IDs (index of gene_df).

    Returns (symbol -> ensembl_id, list of missing symbols).
    """
    if symbol_col not in gene_df.columns:
        raise KeyError(f"gene metadata missing column {symbol_col!r}")

    symbol_to_ensembl: dict[str, str] = {}
    missing: list[str] = []

    for sym in symbols:
        mask = gene_df[symbol_col] == sym
        hits = gene_df.index[mask].tolist()
        if not hits:
            missing.append(sym)
        elif len(hits) > 1:
            warnings.warn(
                f"Gene symbol {sym!r} maps to {len(hits)} Ensembl IDs; using first.",
                UserWarning,
                stacklevel=2,
            )
            symbol_to_ensembl[sym] = hits[0]
        else:
            symbol_to_ensembl[sym] = hits[0]

    return symbol_to_ensembl, missing


def warn_missing_genes(found: list[str], requested: list[str]) -> None:
    """Warn for genes not found in the dataset; do not raise."""
    missing = set(requested) - set(found)
    if missing:
        warnings.warn(
            f"Genes not found in dataset (skipped): {sorted(missing)}",
            UserWarning,
            stacklevel=2,
        )


def scrna_pooled_column_label(dissection: str, requested_areas: list[str]) -> str:
    """Heatmap column label for a dissection ROI pooling multiple config areas."""
    return f"{dissection} (pooled: {', '.join(requested_areas)})"


def scrna_heatmap_columns(config: dict[str, Any]) -> list[str]:
    """Display column labels for scRNA heatmaps (one column per dissection pool)."""
    brain_areas: list[str] = config["brain_areas"]
    pools: dict[str, list[str]] = config.get("_scrna_pools") or {}
    pooled_areas = {a for areas in pools.values() for a in areas}

    columns: list[str] = []
    seen_dissections: set[str] = set()
    for ba in brain_areas:
        if ba in pooled_areas:
            dissection = next(d for d, areas in pools.items() if ba in areas)
            if dissection in seen_dissections:
                continue
            seen_dissections.add(dissection)
            columns.append(scrna_pooled_column_label(dissection, pools[dissection]))
        else:
            columns.append(ba)
    return columns


def scrna_column_to_brain_area(config: dict[str, Any]) -> dict[str, str]:
    """Map heatmap display column label -> brain_area value in aggregated data."""
    pools: dict[str, list[str]] = config.get("_scrna_pools") or {}
    pooled_areas = {a for areas in pools.values() for a in areas}

    mapping: dict[str, str] = {}
    seen_dissections: set[str] = set()
    for ba in config["brain_areas"]:
        if ba in pooled_areas:
            dissection = next(d for d, areas in pools.items() if ba in areas)
            if dissection in seen_dissections:
                continue
            seen_dissections.add(dissection)
            label = scrna_pooled_column_label(dissection, pools[dissection])
            mapping[label] = dissection
        else:
            mapping[ba] = ba
    return mapping


def build_brain_area_mapping(
    cache: Any,
    brain_areas: list[str],
) -> tuple[
    dict[str, set[str]],
    Callable[[pd.DataFrame], pd.Series],
    dict[str, list[str]],
]:
    """
    Build region_of_interest_acronym -> brain_area mapping for WMB-10X scRNA.

    When multiple config brain_areas share one dissection ROI (e.g. VISp/VISpm/VISam
    -> VIS), cells are labeled with the dissection acronym and the returned pool map
    records which config areas were collapsed.

    Returns (area_to_rois, assign_brain_area, scrna_pools).
    """
    roi_meta = cache.get_metadata_dataframe(
        directory="WMB-10X",
        file_name="region_of_interest_metadata",
    )
    all_rois = set(roi_meta["acronym"].astype(str))

    roi_to_area: dict[str, str] = {}
    area_to_rois: dict[str, set[str]] = {ba: set() for ba in brain_areas}
    scrna_pools: dict[str, list[str]] = {}

    dissection_to_config: dict[str, list[str]] = {}
    for ba in brain_areas:
        dissection = _CCF_TO_DISSECTION_ROI.get(ba)
        if dissection and dissection in all_rois:
            dissection_to_config.setdefault(dissection, []).append(ba)

    for dissection, config_areas in dissection_to_config.items():
        if len(config_areas) > 1:
            scrna_pools[dissection] = list(config_areas)
            roi_to_area[dissection] = dissection
            for ba in config_areas:
                area_to_rois[ba].add(dissection)
            warnings.warn(
                f"brain_areas {config_areas} share WMB-10X dissection ROI {dissection!r} "
                f"(scRNA-seq has no CCF parcellation). Cells labeled {dissection!r}; "
                f"heatmap shows one pooled column. Use MERFISH for sub-region resolution.",
                UserWarning,
                stacklevel=2,
            )
        else:
            ba = config_areas[0]
            roi_to_area[dissection] = ba
            area_to_rois[ba].add(dissection)
            warnings.warn(
                f"brain_area {ba!r} maps to WMB-10X dissection ROI {dissection!r} "
                f"(scRNA-seq has no CCF parcellation at this resolution).",
                UserWarning,
                stacklevel=2,
            )

    for roi in all_rois:
        if roi in roi_to_area:
            continue
        area = _roi_to_brain_area(roi)
        if area is not None and area in brain_areas:
            roi_to_area[roi] = area
            area_to_rois[area].add(roi)

    # Also map via feature_matrix_label for cells whose ROI might be ambiguous.
    def assign_brain_area(cell_meta: pd.DataFrame) -> pd.Series:
        roi_col = "region_of_interest_acronym"
        pkg_col = "feature_matrix_label"

        def _assign_row(row: pd.Series) -> str | float:
            roi = str(row.get(roi_col, ""))
            if roi in roi_to_area:
                return roi_to_area[roi]
            pkg = str(row.get(pkg_col, ""))
            if pkg.startswith("WMB-10Xv3-"):
                suffix = pkg.replace("WMB-10Xv3-", "")
                area = _PACKAGE_TO_BRAIN_AREA.get(suffix)
                if area and area in brain_areas:
                    return area
            return np.nan

        return cell_meta.apply(_assign_row, axis=1)

    return area_to_rois, assign_brain_area, scrna_pools


def _roi_to_brain_area(roi: str) -> str | None:
    """Map a single dissection ROI acronym to a config brain_area."""
    if roi in ("STRd", "STRv"):
        return "STR"
    if roi == "sAMY":
        return "AMY"
    if roi in _CORTICAL_ROIS:
        return "CTX"
    # Direct matches (same acronym as CCF division in config).
    direct = {"TH", "HY", "MB", "MY", "CB", "PAL", "P", "HIP", "CTX", "STR", "AMY"}
    if roi in direct:
        return roi
    if roi == "ENT":
        return "HIP"
    if roi == "RHP":
        return "HIP"
    if roi == "LSX":
        return "PAL"
    if roi == "CTXsp":
        return "CTX"
    return None


def filter_cell_types_by_name(
    cell_types: list[str] | pd.Index,
    config: dict[str, Any],
) -> list[str]:
    """
    Apply optional ``cell_type_name_filter`` from config (substring match, OR logic).

    Returns sorted cell type names; empty filter keeps all.
    """
    patterns: list[str] = config.get("cell_type_name_filter") or []
    names = [str(ct) for ct in cell_types]
    if not patterns:
        return sorted(names)
    matched = sorted(ct for ct in names if any(p in ct for p in patterns))
    if not matched:
        warnings.warn(
            f"cell_type_name_filter {patterns!r} matched no cell types "
            f"(of {len(names)} at {config.get('cell_type_level', 'unknown')!r})",
            UserWarning,
            stacklevel=2,
        )
    return matched


def print_path(prefix: str, path: str | Path) -> None:
    """Print a filesystem path; in Jupyter, render as a clickable ``file://`` link.

    Plain ``print`` of paths containing ``!`` breaks notebook auto-linking (e.g.
    ``…/!Projects/…``); explicit HTML links avoid that.
    """
    p = Path(path).expanduser().resolve()
    uri = p.as_uri()
    try:
        from IPython import get_ipython
        from IPython.display import HTML, display

        if get_ipython() is not None:
            display(HTML(f'{prefix} <a href="{uri}">{p}</a>'))
            return
    except ImportError:
        pass
    print(f"{prefix} {p}")


def top_variable_cell_types(
    matrix: pd.DataFrame,
    n: int = 50,
) -> list[str]:
    """
    Return top-N cell types (index) by variance across columns (genes/regions).

    matrix: rows = cell types, columns = genes or region columns.
    """
    if matrix.empty or len(matrix.index) < 2:
        return list(matrix.index)
    variances = matrix.var(axis=1)
    variances = variances.dropna()
    if variances.empty:
        return list(matrix.index[:n])
    return variances.nlargest(min(n, len(variances))).index.tolist()
