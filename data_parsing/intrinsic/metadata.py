"""Metadata joins from All cells and metadata-only sheets."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from openpyxl import Workbook

from .area import resolve_area
from .config import ALL_CELLS_SHEET, CLUSTER_SHEET, METADATA_SHEETS
from .ids import normalize_id


@dataclass
class CellMetadata:
    source_sheet: str = ""
    region: str = ""
    layer: str = ""
    classic_burster: Any = None
    area: str = ""
    area_ccf: str = ""
    area_mismatch: bool = False
    exclude_flag: Any = None
    cre_label: Any = None
    axon: Any = None
    note: Any = None
    comment: Any = None
    excluded_in_may: Any = None
    time_from_5ht: Any = None
    layer_meta: Any = None
    assumed_type: str = ""
    projection_target: str = ""
    physiological_cluster: str = ""
    task_note: str = ""
    caesum_note: str = ""
    dup_conflict: bool = False


def _normalize_pt(val: str) -> str:
    return "ET" if val.strip().upper() == "PT" else val


def load_all_cells(wb: Workbook) -> dict[str, dict[str, Any]]:
    ws = wb[ALL_CELLS_SHEET]
    out: dict[str, dict[str, Any]] = {}
    for r in range(1, ws.max_row + 1):
        cid_raw = ws.cell(r, 1).value
        if cid_raw is None:
            continue
        cid = normalize_id(str(cid_raw).strip())
        if not cid.startswith("nm"):
            continue
        out[cid] = {
            "exclude_flag": ws.cell(r, 2).value,
            "time_from_5ht": ws.cell(r, 4).value,
            "cre_label": ws.cell(r, 7).value,
            "all_cells_area": ws.cell(r, 8).value,
            "axon": ws.cell(r, 9).value,
            "note": ws.cell(r, 10).value,
            "layer_meta": ws.cell(r, 11).value,
            "excluded_in_may": ws.cell(r, 14).value,
            "comment": ws.cell(r, 15).value,
        }
    return out


def load_metadata_tags(wb: Workbook) -> dict[str, dict[str, str]]:
    """IDs → tag dict from metadata-only sheets."""
    tags: dict[str, dict[str, str]] = {}
    for sheet_name, tag_values in METADATA_SHEETS.items():
        if sheet_name not in wb.sheetnames:
            continue
        ws = wb[sheet_name]
        for c in range(1, ws.max_column + 1):
            v = ws.cell(1, c).value
            if v is None:
                continue
            s = str(v).strip()
            if not (s.lower().startswith("nm") or __import__("re").match(r"^\d{4}_", s)):
                continue
            cid = normalize_id(s)
            tags.setdefault(cid, {}).update(tag_values)
    return tags


def load_cluster_metadata(wb: Workbook) -> dict[str, str]:
    if CLUSTER_SHEET not in wb.sheetnames:
        return {}
    ws = wb[CLUSTER_SHEET]
    out: dict[str, str] = {}
    for r in range(1, ws.max_row + 1):
        cluster = ws.cell(r, 9).value  # col I
        cid_raw = ws.cell(r, 10).value  # col J
        if cid_raw is None or cluster is None:
            continue
        cid = normalize_id(str(cid_raw).strip())
        out[cid] = str(cluster).strip()
    return out


def fill_assumed_type_from_cluster(meta: CellMetadata) -> bool:
    """Fill assumed_type from physiological_cluster when empty. Returns True if filled."""
    if meta.assumed_type or not meta.physiological_cluster:
        return False
    cl = meta.physiological_cluster
    if cl == "IT":
        meta.assumed_type = "IT"
    elif cl in ("ET1", "ET2"):
        meta.assumed_type = "ET"
    return True


def merge_cell_metadata(
    cell_id: str,
    *,
    source_sheet: str,
    region: str,
    layer: str,
    sheet_meta: dict[str, Any],
    all_cells: dict[str, dict[str, Any]],
    tags: dict[str, dict[str, str]],
    clusters: dict[str, str],
    task_note: str = "",
    caesum_note: str = "",
) -> CellMetadata:
    ac = all_cells.get(cell_id, {})
    morph_raw = sheet_meta.get("area_morph_raw", {}).get(cell_id)
    area, area_ccf, mismatch = resolve_area(
        sheet_region=region,
        morph_raw=morph_raw,
        all_cells_area=ac.get("all_cells_area"),
        all_cells_note=ac.get("note"),
    )

    meta = CellMetadata(
        source_sheet=source_sheet,
        region=region,
        layer=layer,
        classic_burster=sheet_meta.get("classic_burster", {}).get(cell_id),
        area=area,
        area_ccf=area_ccf,
        area_mismatch=mismatch,
        exclude_flag=ac.get("exclude_flag"),
        cre_label=ac.get("cre_label"),
        axon=ac.get("axon"),
        note=ac.get("note"),
        comment=ac.get("comment"),
        excluded_in_may=ac.get("excluded_in_may"),
        time_from_5ht=ac.get("time_from_5ht"),
        layer_meta=ac.get("layer_meta"),
        task_note=task_note,
        caesum_note=caesum_note,
    )

    for k, v in tags.get(cell_id, {}).items():
        setattr(meta, k if k != "assumed_type" else "assumed_type", v)

    if cell_id in clusters:
        meta.physiological_cluster = clusters[cell_id]

    if meta.assumed_type:
        meta.assumed_type = _normalize_pt(meta.assumed_type)

    fill_assumed_type_from_cluster(meta)
    return meta


def metadata_to_row(meta: CellMetadata) -> dict[str, Any]:
    return {
        "source_sheet": meta.source_sheet,
        "region": meta.region,
        "layer": meta.layer,
        "classic_burster": meta.classic_burster,
        "area": meta.area,
        "areaCCF": meta.area_ccf,
        "area_mismatch": meta.area_mismatch,
        "exclude_flag": meta.exclude_flag,
        "cre_label": meta.cre_label,
        "axon": meta.axon,
        "note": meta.note,
        "comment": meta.comment,
        "excluded_in_May": meta.excluded_in_may,
        "time_from_5HT": meta.time_from_5ht,
        "layer_meta": meta.layer_meta,
        "assumed_type": meta.assumed_type,
        "projection_target": meta.projection_target,
        "physiological_cluster": meta.physiological_cluster,
        "task_note": meta.task_note,
        "caesum_note": meta.caesum_note,
        "dup_conflict": meta.dup_conflict,
    }
