"""Gene resolution, brain-area mapping, and helper utilities."""

from __future__ import annotations

import warnings
from typing import Any, Callable

import numpy as np
import pandas as pd

# Cortical dissection ROI acronyms (WMB-10X region_of_interest_acronym).
_CORTICAL_ROIS = frozenset({
    "ACA", "AI", "AUD", "AUD-TEa-PERI-ECT", "MO-FRP", "MOp", "PL-ILA-ORB",
    "RSP", "SS-GU-VISC", "SSp", "TEa-PERI-ECT", "VIS", "VIS-PTLp",
})

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


def build_brain_area_mapping(
    cache: Any,
    brain_areas: list[str],
) -> tuple[dict[str, set[str]], Callable[[pd.DataFrame], pd.Series]]:
    """
    Build region_of_interest_acronym -> brain_area mapping for WMB-10X scRNA.

    Returns (brain_area -> set of ROI acronyms), assign_brain_area function).
    """
    roi_meta = cache.get_metadata_dataframe(
        directory="WMB-10X",
        file_name="region_of_interest_metadata",
    )
    all_rois = set(roi_meta["acronym"].astype(str))

    roi_to_area: dict[str, str] = {}
    area_to_rois: dict[str, set[str]] = {ba: set() for ba in brain_areas}

    for roi in all_rois:
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

    return area_to_rois, assign_brain_area


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
