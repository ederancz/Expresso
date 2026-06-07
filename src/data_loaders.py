"""ABC Atlas data loading with backed h5ad partial reads."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from tqdm import tqdm

from src.config import get_cache_dir, get_expression_suffix
from src.utils import (
    build_brain_area_mapping,
    resolve_gene_ids,
    warn_missing_genes,
)

WMB_10X_DIR = "WMB-10X"
WMB_10Xv3_DIR = "WMB-10Xv3"
WMB_TAXONOMY_DIR = "WMB-taxonomy"
MERFISH_DIR = "MERFISH-C57BL6J-638850"
MERFISH_IMPUTED_DIR = "MERFISH-C57BL6J-638850-imputed"
MERFISH_CCF_DIR = "MERFISH-C57BL6J-638850-CCF"


def get_abc_cache(config: dict[str, Any]) -> Any:
    """Initialise AbcProjectCache from config cache_dir."""
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

    cache_dir = get_cache_dir(config)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return AbcProjectCache.from_cache_dir(cache_dir)


def load_scrna_cell_metadata(cache: Any, config: dict[str, Any]) -> pd.DataFrame:
    """
    Load WMB-10Xv3 cells with taxonomy and brain_area assignment.

    Filtered to config brain_areas. Indexed by cell_label.
    """
    cell = cache.get_metadata_dataframe(
        directory=WMB_10X_DIR,
        file_name="cell_metadata",
        dtype={"cell_label": str},
    )
    cell = cell.set_index("cell_label")
    cell = cell[cell["dataset_label"] == "WMB-10Xv3"].copy()

    cluster_details = cache.get_metadata_dataframe(
        directory=WMB_TAXONOMY_DIR,
        file_name="cluster_to_cluster_annotation_membership_pivoted",
        keep_default_na=False,
    )
    cluster_details = cluster_details.set_index("cluster_alias")

    cell_type_col = config["cell_type_level"]
    brain_areas = config["brain_areas"]

    area_to_rois, assign_brain_area, scrna_pools = build_brain_area_mapping(
        cache, brain_areas,
    )
    config["_scrna_pools"] = scrna_pools

    cell = cell.join(cluster_details, on="cluster_alias")
    cell["brain_area"] = assign_brain_area(cell)

    assignable = set(brain_areas) | set(scrna_pools.keys())
    cell = cell[cell["brain_area"].isin(assignable)].copy()

    pooled_config_areas = {a for areas in scrna_pools.values() for a in areas}
    for ba in brain_areas:
        if ba in pooled_config_areas:
            continue
        n = (cell["brain_area"] == ba).sum()
        if n == 0:
            warnings.warn(
                f"No cells mapped to brain_area {ba!r}; "
                f"ROI set was {sorted(area_to_rois.get(ba, []))}",
                UserWarning,
                stacklevel=2,
            )

    for dissection, config_areas in scrna_pools.items():
        n = (cell["brain_area"] == dissection).sum()
        if n == 0:
            warnings.warn(
                f"No cells mapped to pooled dissection ROI {dissection!r} "
                f"(config areas {config_areas}); "
                f"ROI set was {sorted(area_to_rois.get(config_areas[0], []))}",
                UserWarning,
                stacklevel=2,
            )

    if cell_type_col not in cell.columns:
        raise KeyError(f"cell_type_level {cell_type_col!r} not in joined metadata")

    return cell


def check_scrna_genes_in_metadata(
    cache: Any,
    genes: list[str],
) -> tuple[list[str], list[str]]:
    """
    Check gene symbols against WMB-10X gene metadata.

    Returns (found_symbols, missing_symbols).
    """
    gene_df = cache.get_metadata_dataframe(directory=WMB_10X_DIR, file_name="gene")
    gene_df = gene_df.set_index("gene_identifier")
    symbol_to_ensembl, missing = resolve_gene_ids(gene_df, genes)
    return list(symbol_to_ensembl.keys()), missing


def load_expression_subset(
    cache: Any,
    genes: list[str],
    cell_meta: pd.DataFrame,
    config: dict[str, Any],
) -> ad.AnnData | None:
    """
    Load expression for genes x filtered cells without reading full matrices.

    Iterates WMB-10Xv3 packages (feature_matrix_label) present in cell_meta.
    """
    if not genes or cell_meta.empty:
        return None

    gene_df = cache.get_metadata_dataframe(directory=WMB_10X_DIR, file_name="gene")
    gene_df = gene_df.set_index("gene_identifier")

    symbol_to_ensembl, missing = resolve_gene_ids(gene_df, genes)
    warn_missing_genes(list(symbol_to_ensembl), genes)

    if not symbol_to_ensembl:
        return None

    ensembl_ids = list(symbol_to_ensembl.values())
    ensembl_to_symbol = {v: k for k, v in symbol_to_ensembl.items()}

    data_suffix = get_expression_suffix(config)
    packages = cell_meta["feature_matrix_label"].dropna().unique()

    frames: list[ad.AnnData] = []
    for pkg in tqdm(packages, desc="WMB-10Xv3 packages"):
        cells_in_pkg = cell_meta.index[cell_meta["feature_matrix_label"] == pkg]
        if len(cells_in_pkg) == 0:
            continue

        file_path = cache.get_file_path(
            directory=WMB_10Xv3_DIR,
            file_name=f"{pkg}/{data_suffix}",
        )
        adata = ad.read_h5ad(file_path, backed="r")
        try:
            cells_present = adata.obs_names.intersection(cells_in_pkg)
            if len(cells_present) == 0:
                continue

            genes_in_file = [g for g in ensembl_ids if g in adata.var_names]
            if not genes_in_file:
                continue

            subset = adata[cells_present, genes_in_file].to_memory()
            frames.append(subset)
        finally:
            adata.file.close()

    if not frames:
        return None

    combined = ad.concat(frames, join="outer")
    combined.var["gene_symbol"] = [
        ensembl_to_symbol.get(e, gene_df.loc[e, "gene_symbol"] if e in gene_df.index else e)
        for e in combined.var_names
    ]
    return combined


def aggregate_scrna_expression(
    adata: ad.AnnData,
    cell_meta: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Mean expression per (cell_type_level, brain_area).

    Returns long DataFrame: cell_type, brain_area, gene, mean_expression, family.
    """
    cell_type_col = config["cell_type_level"]
    genes_flat = config["_genes_flat"]

    if hasattr(adata.X, "toarray"):
        x = adata.X.toarray()
    else:
        x = np.asarray(adata.X)

    gene_symbols = adata.var["gene_symbol"].tolist()
    df = pd.DataFrame(x, index=adata.obs_names, columns=gene_symbols)

    meta = cell_meta.loc[df.index, [cell_type_col, "brain_area"]]
    df = df.join(meta)

    grouped = df.groupby([cell_type_col, "brain_area"], observed=True)[gene_symbols].mean()

    rows = []
    for (ct, ba), series in grouped.iterrows():
        for gene, val in series.items():
            rows.append({
                "cell_type": ct,
                "brain_area": ba,
                "gene": gene,
                "mean_expression": val,
                "family": genes_flat.get(gene, "unknown"),
            })

    return pd.DataFrame(rows)


def combined_heatmap_matrix(
    agg_long: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Combined heatmap matrix: rows=cell types (alphabetical), columns=receptor genes
    (family order).

    Expression is mean across configured brain areas per cell type × gene.
    Cell types may be narrowed via ``cell_type_name_filter`` in config.
    """
    from src.utils import filter_cell_types_by_name

    genes = config["_all_genes"]
    cell_types = filter_cell_types_by_name(agg_long["cell_type"].unique(), config)
    mat = agg_long.pivot_table(
        index="cell_type",
        columns="gene",
        values="mean_expression",
        aggfunc="mean",
    )
    return mat.reindex(index=cell_types, columns=genes)


def family_gene_region_matrix(
    agg_long: pd.DataFrame,
    family: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Per-family heatmap: rows=cell types, columns=brain areas.

    For multiple genes, average expression across genes in the family.
    Pooled scRNA dissection ROIs appear as a single column (see scrna_heatmap_columns).
    """
    from src.utils import (
        filter_cell_types_by_name,
        scrna_column_to_brain_area,
        scrna_heatmap_columns,
    )

    sub = agg_long[agg_long["family"] == family]
    if sub.empty:
        return pd.DataFrame()
    mat = sub.pivot_table(
        index="cell_type",
        columns="brain_area",
        values="mean_expression",
        aggfunc="mean",
    )
    col_map = scrna_column_to_brain_area(config)
    display_cols = scrna_heatmap_columns(config)

    pools: dict[str, list[str]] = config.get("_scrna_pools") or {}

    renamed = pd.DataFrame(index=mat.index)
    for display, source_ba in col_map.items():
        if source_ba in mat.columns:
            renamed[display] = mat[source_ba]
        else:
            # Fallback: stale agg_long may still use pre-pool labels (e.g. VISam).
            legacy_cols = [c for c in pools.get(source_ba, []) if c in mat.columns]
            if legacy_cols:
                renamed[display] = mat[legacy_cols].mean(axis=1)
            else:
                renamed[display] = np.nan
    cell_types = filter_cell_types_by_name(renamed.index, config)
    return renamed.reindex(index=cell_types, columns=display_cols)


def load_merfish_cell_metadata(cache: Any, config: dict[str, Any]) -> pd.DataFrame:
    """MERFISH cells with cluster annotation, section coords, and CCF coords."""
    dataset = config["data"].get("merfish_dataset", MERFISH_DIR)

    cell = cache.get_metadata_dataframe(
        directory=dataset,
        file_name="cell_metadata_with_cluster_annotation",
        dtype={"cell_label": str},
    )
    cell = cell.set_index("cell_label")

    # Section coordinates (retained but not used for M2 plots).
    cell = cell.rename(columns={"x": "x_section", "y": "y_section", "z": "z_section"})

    ccf_coords = cache.get_metadata_dataframe(
        directory=MERFISH_CCF_DIR,
        file_name="ccf_coordinates",
        dtype={"cell_label": str},
    )
    ccf_coords = ccf_coords.set_index("cell_label")
    ccf_coords = ccf_coords.rename(
        columns={"x": "x_ccf", "y": "y_ccf", "z": "z_ccf"},
    )
    if "parcellation_index" in ccf_coords.columns:
        ccf_coords = ccf_coords.drop(columns=["parcellation_index"])

    cell = cell.join(ccf_coords, how="inner")
    return cell


def _merfish_gene_panel(cache: Any, directory: str) -> set[str]:
    gene_df = cache.get_metadata_dataframe(directory=directory, file_name="gene")
    return set(gene_df["gene_symbol"].astype(str))


def _imputed_gene_symbols(cache: Any) -> set[str]:
    """Load gene symbols from imputed h5ad var (metadata only via backed read)."""
    data_suffix = "log2"
    path = cache.get_file_path(
        directory=MERFISH_IMPUTED_DIR,
        file_name=f"C57BL6J-638850-imputed/{data_suffix}",
    )
    adata = ad.read_h5ad(path, backed="r")
    try:
        if "gene_symbol" in adata.var.columns:
            return set(adata.var["gene_symbol"].astype(str))
        return set()
    finally:
        adata.file.close()


_imputed_panel_cache: set[str] | None = None


def check_gene_availability(
    cache: Any,
    gene: str,
    config: dict[str, Any],
) -> str:
    """Return 'present', 'imputed', or 'missing'."""
    global _imputed_panel_cache

    dataset = config["data"].get("merfish_dataset", MERFISH_DIR)
    panel = _merfish_gene_panel(cache, dataset)
    if gene in panel:
        return "present"

    if not config["data"].get("use_imputed_merfish", True):
        return "missing"

    if _imputed_panel_cache is None:
        try:
            _imputed_panel_cache = _imputed_gene_symbols(cache)
        except Exception:
            _imputed_panel_cache = set()

    if gene in _imputed_panel_cache:
        return "imputed"
    return "missing"


def load_single_gene_merfish(
    cache: Any,
    gene: str,
    config: dict[str, Any],
) -> tuple[pd.Series | None, str]:
    """
    Load one gene's expression for all MERFISH cells.

    Returns (Series indexed by cell_label, source_label).
    source_label: 'measured', 'imputed', or 'missing'.
    """
    status = check_gene_availability(cache, gene, config)
    if status == "missing":
        return None, "missing"

    data_suffix = get_expression_suffix(config)
    if status == "present":
        directory = config["data"].get("merfish_dataset", MERFISH_DIR)
        file_name = f"C57BL6J-638850/{data_suffix}"
    else:
        directory = MERFISH_IMPUTED_DIR
        file_name = f"C57BL6J-638850-imputed/{data_suffix}"

    path = cache.get_file_path(directory=directory, file_name=file_name)
    adata = ad.read_h5ad(path, backed="r")
    try:
        if "gene_symbol" not in adata.var.columns:
            return None, "missing"

        mask = adata.var["gene_symbol"] == gene
        if not mask.any():
            return None, "missing"

        ensembl_id = adata.var.index[mask][0]
        col = adata[:, ensembl_id]
        if hasattr(col.X, "toarray"):
            expr = np.asarray(col.X.toarray()).ravel()
        else:
            expr = np.asarray(col.X).ravel()

        source = "measured" if status == "present" else "imputed"
        return pd.Series(expr, index=adata.obs_names, name=gene), source
    finally:
        adata.file.close()
