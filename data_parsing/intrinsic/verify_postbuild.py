"""Phase 3 post-build verification — prove Prime Directive held after assembly."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from .area import resolve_area
from .config import DRUG_SHEETS
from .dedup import merge_sheet_metas, pick_canonical
from .ids import HeaderInfo
from .labels import col_name
from .sheet_parser import SheetParseResult, sheet_region_layer
from .verify import MAX_BYTE_EQUAL_ERRORS, byte_equal

CESIUM_SHEET = "V2M_L5_Caesum"
TASK_SHEET = "V2M_L5_TASK_Acidic_pH"
TASK_EXCLUDED_CELL = "nm2025_06_18_c1"


class PostBuildVerificationError(Exception):
    """Raised when Phase 3 post-build checks fail."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        super().__init__(self.summary())

    def summary(self) -> str:
        head = f"Phase 3 post-build verification failed ({len(self.errors)} error(s))"
        if not self.errors:
            return head
        shown = self.errors[:MAX_BYTE_EQUAL_ERRORS]
        tail = ""
        if len(self.errors) > MAX_BYTE_EQUAL_ERRORS:
            tail = f"\n  … and {len(self.errors) - MAX_BYTE_EQUAL_ERRORS} more"
        return head + ":\n  " + "\n  ".join(shown) + tail


def _ordered_param_columns(all_params: set[str]) -> list[str]:
    section_order = ["_IV", "_short_depol", "_EPSP", "_crit_freq", "_sag", "_chirp"]
    by_sec: dict[str, list[str]] = defaultdict(list)
    for p in all_params:
        sec = p.split("__", 1)[0]
        by_sec[sec].append(p)
    out: list[str] = []
    for sec in section_order:
        out.extend(sorted(by_sec.get(sec, [])))
    for sec in sorted(by_sec):
        if sec not in section_order:
            out.extend(sorted(by_sec[sec]))
    return out


def _expected_param_columns(parse_results: list[SheetParseResult]) -> list[str]:
    cols: set[str] = set()
    for pr in parse_results:
        for pv in pr.values:
            cols.add(col_name(pv.section, pv.label))
    return _ordered_param_columns(cols)


def _merge_notes(note: Any, comment: Any) -> str:
    parts: list[str] = []
    for val in (note, comment):
        if val is None:
            continue
        text = str(val).strip()
        if text:
            parts.append(text)
    return " | ".join(parts)


def _is_excluded(cid: str, pr: SheetParseResult, dropped_ids: set[str]) -> bool:
    return cid in pr.task_exclude or cid in dropped_ids


def _expected_output_ids(
    parse_results: list[SheetParseResult],
    dropped_ids: set[str],
) -> tuple[set[str], set[tuple[str, str]]]:
    """Return (expected control cell_ids, expected effect (cell_id, experiment) pairs)."""
    control: set[str] = set()
    effect: set[tuple[str, str]] = set()

    for pr in parse_results:
        sheet = pr.sheet_name
        cells_with_control: set[str] = set()
        cells_with_effect: set[str] = set()
        for pv in pr.values:
            if _is_excluded(pv.cell_id, pr, dropped_ids):
                continue
            if pv.block == "control":
                cells_with_control.add(pv.cell_id)
            elif pv.block == "effect":
                cells_with_effect.add(pv.cell_id)

        if sheet == CESIUM_SHEET:
            experiment = DRUG_SHEETS[sheet]
            for cid in cells_with_effect:
                effect.add((cid, experiment))
        elif sheet in DRUG_SHEETS:
            experiment = DRUG_SHEETS[sheet]
            control |= cells_with_control
            for cid in cells_with_effect:
                effect.add((cid, experiment))
        else:
            control |= cells_with_control

    return control, effect


def _all_header_ids(parse_results: list[SheetParseResult]) -> dict[str, list[tuple[str, str]]]:
    """Map normalized header ID → list of (sheet_name, raw_header)."""
    out: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for pr in parse_results:
        for h in pr.headers:
            out[h.normalized].append((pr.sheet_name, h.raw))
    return out


def verify_cell_id_fidelity(
    parse_results: list[SheetParseResult],
    *,
    control_rows: list[dict[str, Any]],
    effect_rows: list[dict[str, Any]],
    dropped_ids: set[str],
) -> list[str]:
    errors: list[str] = []
    expected_control, expected_effect = _expected_output_ids(parse_results, dropped_ids)

    actual_control = {r["cell_id"] for r in control_rows}
    actual_effect = {(r["cell_id"], r["experiment"]) for r in effect_rows}

    for cid in sorted(expected_control - actual_control):
        errors.append(f"cell-id: missing from control_excitability: {cid!r}")
    for cid in sorted(actual_control - expected_control):
        errors.append(f"cell-id: unexpected in control_excitability: {cid!r}")

    for pair in sorted(expected_effect - actual_effect):
        errors.append(
            f"cell-id: missing from pharmacology_effect: {pair[0]!r} × {pair[1]!r}"
        )
    for pair in sorted(actual_effect - expected_effect):
        errors.append(
            f"cell-id: unexpected in pharmacology_effect: {pair[0]!r} × {pair[1]!r}"
        )

    output_ids = actual_control | {cid for cid, _ in actual_effect}
    if TASK_EXCLUDED_CELL in output_ids:
        errors.append(
            f"cell-id: TASK-excluded {TASK_EXCLUDED_CELL!r} must not appear in output"
        )

    header_ids = _all_header_ids(parse_results)
    for cid in sorted(output_ids):
        lookup = cid.split("#", 1)[0] if "#" in cid else cid
        if lookup not in header_ids and cid not in header_ids:
            errors.append(f"cell-id: output ID {cid!r} has no source header")

    for pr in parse_results:
        for h in pr.headers:
            if not h.normalized.startswith("nm"):
                errors.append(
                    f"cell-id: nm-prefix: {pr.sheet_name} col {h.col} "
                    f"raw={h.raw!r} → {h.normalized!r}"
                )
            if _is_excluded(h.normalized, pr, dropped_ids):
                if h.normalized in actual_control:
                    errors.append(
                        f"cell-id: excluded {h.normalized!r} on {pr.sheet_name} "
                        f"present in control_excitability"
                    )
                if any(cid == h.normalized for cid, _ in actual_effect):
                    errors.append(
                        f"cell-id: excluded {h.normalized!r} on {pr.sheet_name} "
                        f"present in pharmacology_effect"
                    )

        if pr.sheet_name == TASK_SHEET:
            raw_counts: dict[str, list[HeaderInfo]] = defaultdict(list)
            for h in pr.headers:
                base = h.normalized.split("#")[0]
                raw_counts[base].append(h)
            for base, headers in raw_counts.items():
                if len(headers) > 1:
                    suffixes = sorted(
                        h.normalized.split("#", 1)[1]
                        for h in headers
                        if "#" in h.normalized
                    )
                    if suffixes != ["1", "2"]:
                        errors.append(
                            f"cell-id: TASK duplicate {base!r} must disambiguate as "
                            f"#1/#2; got {[h.normalized for h in headers]!r}"
                        )
                    for h in headers:
                        if not h.task_note:
                            errors.append(
                                f"cell-id: TASK header {h.normalized!r} missing task_note"
                            )

    return errors



def verify_numerical_fidelity(
    parse_results: list[SheetParseResult],
    *,
    control_rows: list[dict[str, Any]],
    effect_rows: list[dict[str, Any]],
    dropped_ids: set[str],
) -> tuple[list[str], int]:
    errors: list[str] = []
    control_by_id = {r["cell_id"]: r for r in control_rows}
    effect_by_key = {(r["cell_id"], r["experiment"]): r for r in effect_rows}
    checked = 0

    for pr in parse_results:
        sheet = pr.sheet_name
        for pv in pr.values:
            if _is_excluded(pv.cell_id, pr, dropped_ids):
                continue
            if pv.block == "control" and sheet == CESIUM_SHEET:
                continue

            row: dict[str, Any] | None = None
            if pv.block == "control":
                row = control_by_id.get(pv.cell_id)
            elif pv.block == "effect" and sheet in DRUG_SHEETS:
                row = effect_by_key.get((pv.cell_id, DRUG_SHEETS[sheet]))

            if row is None:
                errors.append(
                    f"numerical: no output row for {pv.source_sheet} "
                    f"{pv.block} {pv.cell_id!r} "
                    f"{col_name(pv.section, pv.label)!r}"
                )
                continue

            if row.get("source_sheet") != pv.source_sheet:
                continue

            key = col_name(pv.section, pv.label)
            built = row.get(key)
            checked += 1
            if not byte_equal(pv.value, built):
                errors.append(
                    f"numerical: {pv.source_sheet} {pv.block} {pv.cell_id!r} "
                    f"{key!r}: parsed={pv.value!r} built={built!r}"
                )

    return errors, checked


def verify_label_fidelity(
    param_cols: list[str],
    parse_results: list[SheetParseResult],
) -> list[str]:
    expected = _expected_param_columns(parse_results)
    errors: list[str] = []
    if param_cols != expected:
        missing = [c for c in expected if c not in param_cols]
        extra = [c for c in param_cols if c not in expected]
        if missing:
            errors.append(f"label: missing output columns ({len(missing)}): {missing[:5]!r}…")
        if extra:
            errors.append(f"label: extra output columns ({len(extra)}): {extra[:5]!r}…")
        if not missing and not extra:
            errors.append(
                f"label: column order mismatch ({len(param_cols)} cols); "
                f"first diff at index "
                f"{next(i for i, (a, b) in enumerate(zip(param_cols, expected)) if a != b)}"
            )
    return errors


def _parse_instances_for_cell(
    parse_results: list[SheetParseResult],
    cid: str,
) -> list[dict[str, Any]]:
    """Mirror build-time instance list for metadata merge (all sheets carrying this cell)."""
    instances: list[dict[str, Any]] = []
    for pr in parse_results:
        present = any(pv.cell_id == cid for pv in pr.values) or cid in pr.meta.classic_burster
        if not present:
            continue
        region, layer = sheet_region_layer(pr.sheet_name)
        instances.append(
            {
                "source_sheet": pr.sheet_name,
                "region": region,
                "layer": layer,
                "sheet_meta": {
                    "classic_burster": pr.meta.classic_burster,
                    "area_morph_raw": pr.meta.area_morph_raw,
                },
            }
        )
    return instances


def verify_metadata_fidelity(
    parse_results: list[SheetParseResult],
    *,
    control_rows: list[dict[str, Any]],
    effect_rows: list[dict[str, Any]],
    all_cells: dict[str, dict[str, Any]],
    clusters: dict[str, str],
) -> tuple[list[str], int]:
    errors: list[str] = []
    fields_checked = 0

    seen_cluster: set[str] = set()
    for row in control_rows + effect_rows:
        cid = row["cell_id"]
        instances = _parse_instances_for_cell(parse_results, cid)
        merged_meta = merge_sheet_metas(instances) if instances else {}
        canonical = pick_canonical(instances) if instances else {}
        source_sheet = row["source_sheet"]

        expected_cb = merged_meta.get("classic_burster", {}).get(cid)
        fields_checked += 1
        if not byte_equal(row.get("classic_burster"), expected_cb):
            errors.append(
                f"metadata: classic_burster {cid!r} (canonical {source_sheet!r}): "
                f"output={row.get('classic_burster')!r} merged={expected_cb!r}"
            )

        ac = all_cells.get(cid, {})
        morph_raw = merged_meta.get("area_morph_raw", {}).get(cid)
        _area, expected_area_ccf, _mismatch = resolve_area(
            sheet_region=canonical.get("region", ""),
            morph_raw=morph_raw,
            all_cells_area=ac.get("all_cells_area"),
            all_cells_note=ac.get("note"),
        )
        fields_checked += 1
        if row.get("areaCCF") != expected_area_ccf:
            errors.append(
                f"metadata: areaCCF {cid!r} (canonical {source_sheet!r}): "
                f"output={row.get('areaCCF')!r} expected={expected_area_ccf!r} "
                f"(morph={morph_raw!r})"
            )

        for field, ac_key in (
            ("exclude_flag", "exclude_flag"),
            ("cre_label", "cre_label"),
            ("axon", "axon"),
            ("time_from_5HT", "time_from_5ht"),
        ):
            expected = ac.get(ac_key)
            fields_checked += 1
            if not byte_equal(row.get(field), expected):
                errors.append(
                    f"metadata: {field} {cid!r}: "
                    f"output={row.get(field)!r} all_cells={expected!r}"
                )

        expected_notes = _merge_notes(ac.get("note"), ac.get("comment"))
        fields_checked += 1
        if row.get("notes") != expected_notes:
            errors.append(
                f"metadata: notes {cid!r}: "
                f"output={row.get('notes')!r} expected={expected_notes!r}"
            )

        if cid in clusters:
            seen_cluster.add(cid)
            fields_checked += 1
            if row.get("physiological_cluster") != clusters[cid]:
                errors.append(
                    f"metadata: physiological_cluster {cid!r}: "
                    f"output={row.get('physiological_cluster')!r} "
                    f"cluster_analysis_res={clusters[cid]!r}"
                )

    for cid in sorted(clusters):
        if cid in seen_cluster:
            continue
        in_output = any(r["cell_id"] == cid for r in control_rows + effect_rows)
        if in_output:
            errors.append(
                f"metadata: cluster-listed {cid!r} in output but "
                f"physiological_cluster not verified"
            )

    return errors, fields_checked


def run_phase3_checks(
    *,
    parse_results: list[SheetParseResult],
    control_rows: list[dict[str, Any]],
    effect_rows: list[dict[str, Any]],
    param_cols: list[str],
    dropped_ids: set[str],
    all_cells: dict[str, dict[str, Any]],
    clusters: dict[str, str],
) -> dict[str, Any]:
    """
    Run all Phase 3 post-build checks on raw built rows (before sig-fig formatting).
    Raises PostBuildVerificationError on failure.
    Returns a report dict for run_manifest.json.
    """
    cell_errors = verify_cell_id_fidelity(
        parse_results,
        control_rows=control_rows,
        effect_rows=effect_rows,
        dropped_ids=dropped_ids,
    )
    num_errors, values_checked = verify_numerical_fidelity(
        parse_results,
        control_rows=control_rows,
        effect_rows=effect_rows,
        dropped_ids=dropped_ids,
    )
    label_errors = verify_label_fidelity(param_cols, parse_results)
    meta_errors, metadata_fields_checked = verify_metadata_fidelity(
        parse_results,
        control_rows=control_rows,
        effect_rows=effect_rows,
        all_cells=all_cells,
        clusters=clusters,
    )

    errors = cell_errors + num_errors + label_errors + meta_errors
    expected_control, expected_effect = _expected_output_ids(parse_results, dropped_ids)

    report: dict[str, Any] = {
        "passed": len(errors) == 0,
        "cell_id_errors": len(cell_errors),
        "numerical_errors": len(num_errors),
        "label_errors": len(label_errors),
        "metadata_errors": len(meta_errors),
        "values_checked": values_checked,
        "metadata_fields_checked": metadata_fields_checked,
        "expected_control_neurons": len(expected_control),
        "expected_effect_rows": len(expected_effect),
        "actual_control_neurons": len(control_rows),
        "actual_effect_rows": len(effect_rows),
        "param_columns": len(param_cols),
        "task_excluded_cell": TASK_EXCLUDED_CELL,
        "errors": errors,
    }

    if errors:
        raise PostBuildVerificationError(errors)

    return report
