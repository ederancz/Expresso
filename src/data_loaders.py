"""ABC Atlas data loading with backed h5ad partial reads."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from tqdm import tqdm

from scipy.spatial import cKDTree

from src.config import (
    discover_vizgen_samples,
    find_prior_run_parquet,
    get_cache_dir,
    get_expression_suffix,
    get_vizgen_data_dir,
    get_vizgen_samples,
    get_zhuang_datasets,
    vizgen_sample_file_paths,
)
from src.utils import (
    assign_merfish_brain_area,
    build_brain_area_mapping,
    filter_cell_types_by_name,
    resolve_gene_ids,
    warn_missing_genes,
)

WMB_10X_DIR = "WMB-10X"
WMB_10Xv3_DIR = "WMB-10Xv3"
WMB_TAXONOMY_DIR = "WMB-taxonomy"
MERFISH_DIR = "MERFISH-C57BL6J-638850"
MERFISH_IMPUTED_DIR = "MERFISH-C57BL6J-638850-imputed"
MERFISH_CCF_DIR = "MERFISH-C57BL6J-638850-CCF"
ALLEN_CCF_DIR = "Allen-CCF-2020"


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
    Mean expression and detection rate per (cell_type_level, brain_area).

    Returns long DataFrame with columns: cell_type, brain_area, gene,
    mean_expression (mean over all cells in the group, zeros included),
    frac_expressing (fraction of cells with expression > 0), n_cells
    (cells in the group), family.
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

    grouped = df.groupby([cell_type_col, "brain_area"], observed=True)
    mean_df = grouped[gene_symbols].mean()
    frac_df = grouped[gene_symbols].agg(lambda s: float((s > 0).mean()))
    n_cells = grouped.size()

    rows = []
    for (ct, ba), mean_series in mean_df.iterrows():
        frac_series = frac_df.loc[(ct, ba)]
        group_n = int(n_cells.loc[(ct, ba)])
        for gene, val in mean_series.items():
            rows.append({
                "cell_type": ct,
                "brain_area": ba,
                "gene": gene,
                "mean_expression": val,
                "frac_expressing": float(frac_series[gene]),
                "n_cells": group_n,
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


def family_gene_region_matrix_merfish(
    agg_long: pd.DataFrame,
    family: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Per-family MERFISH heatmap: rows=cell types, columns=config brain areas (CCF).

    Unlike scRNA, each config brain area is a separate column (no dissection pooling).
    """
    sub = agg_long[agg_long["family"] == family]
    if sub.empty:
        return pd.DataFrame()
    mat = sub.pivot_table(
        index="cell_type",
        columns="brain_area",
        values="mean_expression",
        aggfunc="mean",
    )
    cell_types = filter_cell_types_by_name(mat.index, config)
    return mat.reindex(index=cell_types, columns=config["brain_areas"])


def _join_merfish_parcellation(
    cache: Any,
    cell: pd.DataFrame,
    dataset: str,
) -> pd.DataFrame:
    """Attach CCF parcellation columns to MERFISH cell metadata."""
    parc_cols = [c for c in cell.columns if c.startswith("parcellation_")]
    if parc_cols:
        return cell

    try:
        parc_cell = cache.get_metadata_dataframe(
            directory=dataset,
            file_name="cell_metadata_with_parcellation_annotation",
            dtype={"cell_label": str},
        )
        parc_cell = parc_cell.set_index("cell_label")
        parc_cols = [
            c for c in parc_cell.columns if c.startswith("parcellation_")
        ]
        if parc_cols:
            return cell.join(parc_cell[parc_cols], how="left")
    except Exception:
        pass

    ccf_coords = cache.get_metadata_dataframe(
        directory=MERFISH_CCF_DIR,
        file_name="ccf_coordinates",
        dtype={"cell_label": str},
    )
    ccf_coords = ccf_coords.set_index("cell_label")
    if "parcellation_index" not in ccf_coords.columns:
        warnings.warn(
            "MERFISH CCF coordinates lack parcellation_index; brain areas unavailable.",
            UserWarning,
            stacklevel=2,
        )
        return cell

    parc_ann = cache.get_metadata_dataframe(
        directory=ALLEN_CCF_DIR,
        file_name="parcellation_to_parcellation_term_membership_acronym",
    )
    parc_ann = parc_ann.set_index("parcellation_index")
    parc_ann.columns = [f"parcellation_{c}" for c in parc_ann.columns]

    cell = cell.join(
        ccf_coords[["parcellation_index"]],
        how="left",
    )
    return cell.join(parc_ann, on="parcellation_index", how="left")


def load_merfish_cell_metadata(cache: Any, config: dict[str, Any]) -> pd.DataFrame:
    """
    MERFISH cells with taxonomy and CCF brain_area assignment.

    Filtered to config brain_areas. Indexed by cell_label.
    """
    dataset = config["data"].get("merfish_dataset", MERFISH_DIR)
    brain_areas = config["brain_areas"]
    cell_type_col = config["cell_type_level"]

    cell = cache.get_metadata_dataframe(
        directory=dataset,
        file_name="cell_metadata_with_cluster_annotation",
        dtype={"cell_label": str},
    )
    cell = cell.set_index("cell_label")
    cell = _join_merfish_parcellation(cache, cell, dataset)
    cell["brain_area"] = assign_merfish_brain_area(cell, brain_areas)
    cell = cell[cell["brain_area"].notna() & cell["brain_area"].isin(brain_areas)].copy()

    for ba in brain_areas:
        n = int((cell["brain_area"] == ba).sum())
        if n == 0:
            warnings.warn(
                f"No MERFISH cells mapped to brain_area {ba!r}",
                UserWarning,
                stacklevel=2,
            )

    if cell_type_col not in cell.columns:
        raise KeyError(f"cell_type_level {cell_type_col!r} not in MERFISH metadata")

    config["_scrna_pools"] = {}
    return cell


def _merfish_gene_panel(cache: Any, directory: str) -> set[str]:
    gene_df = cache.get_metadata_dataframe(directory=directory, file_name="gene")
    return set(gene_df["gene_symbol"].astype(str))


def _imputed_gene_symbols(cache: Any, config: dict[str, Any] | None = None) -> set[str]:
    """Load gene symbols from imputed h5ad var (metadata only via backed read)."""
    data_suffix = get_expression_suffix(config) if config else "log2"
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


def _ensure_float32_x(adata: ad.AnnData) -> ad.AnnData:
    """Cast X to float32; scipy.sparse cannot concat float16 matrices."""
    if hasattr(adata.X, "astype"):
        adata.X = adata.X.astype(np.float32)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)
    return adata


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
            _imputed_panel_cache = _imputed_gene_symbols(cache, config)
        except Exception:
            _imputed_panel_cache = set()

    if gene in _imputed_panel_cache:
        return "imputed"
    return "missing"


def check_merfish_genes(
    cache: Any,
    genes: list[str],
    config: dict[str, Any],
) -> dict[str, list[str]]:
    """Classify genes as measured, imputed, or missing in MERFISH."""
    measured: list[str] = []
    imputed: list[str] = []
    missing: list[str] = []
    for gene in genes:
        status = check_gene_availability(cache, gene, config)
        if status == "present":
            measured.append(gene)
        elif status == "imputed":
            imputed.append(gene)
        else:
            missing.append(gene)
    return {"measured": measured, "imputed": imputed, "missing": missing}


def _load_merfish_h5ad_slice(
    cache: Any,
    config: dict[str, Any],
    genes: list[str],
    cell_ids: pd.Index,
    *,
    imputed: bool,
) -> ad.AnnData | None:
    """Load a genes × cells slice from measured or imputed MERFISH matrix."""
    if not genes:
        return None

    data_suffix = get_expression_suffix(config)
    if imputed:
        directory = MERFISH_IMPUTED_DIR
        file_name = f"C57BL6J-638850-imputed/{data_suffix}"
    else:
        directory = config["data"].get("merfish_dataset", MERFISH_DIR)
        file_name = f"C57BL6J-638850/{data_suffix}"

    path = cache.get_file_path(directory=directory, file_name=file_name)
    adata = ad.read_h5ad(path, backed="r")
    try:
        if "gene_symbol" not in adata.var.columns:
            return None

        cells_present = adata.obs_names.intersection(cell_ids)
        if len(cells_present) == 0:
            return None

        gene_mask = adata.var["gene_symbol"].isin(genes)
        ensembl_ids = adata.var.index[gene_mask].tolist()
        if not ensembl_ids:
            return None

        source_label = "imputed" if imputed else "measured"
        print(
            f"MERFISH {source_label}: reading {len(cells_present):,} cells × "
            f"{len(ensembl_ids)} genes from {path.name} "
            f"(row-first slice; imputed matrix is ~47 GB)…",
            flush=True,
        )
        # CSR-backed h5ad: combined [cells, genes].to_memory() can scan the entire
        # on-disk matrix (hours, 0% CPU). Load filtered rows first, then columns.
        row_subset = adata[cells_present, :].to_memory()
        subset = row_subset[:, ensembl_ids].copy()
        del row_subset
        symbol_by_ensembl = dict(
            zip(subset.var_names, subset.var["gene_symbol"].astype(str))
        )
        subset.var["gene_symbol"] = [
            symbol_by_ensembl[e] for e in subset.var_names
        ]
        return _ensure_float32_x(subset)
    finally:
        adata.file.close()


def load_merfish_expression_subset(
    cache: Any,
    genes: list[str],
    cell_meta: pd.DataFrame,
    config: dict[str, Any],
) -> ad.AnnData | None:
    """
    Batch-load MERFISH expression for requested genes × filtered cells.

    Measured genes come from the ~500-gene panel; imputed genes from the
    ~8k matrix when ``use_imputed_merfish`` is enabled.
    """
    if not genes or cell_meta.empty:
        return None

    availability = check_merfish_genes(cache, genes, config)
    warn_missing_genes(
        availability["measured"] + availability["imputed"],
        genes,
    )

    frames: list[ad.AnnData] = []
    if availability["measured"]:
        print(
            f"Loading {len(availability['measured'])} measured MERFISH genes…",
            flush=True,
        )
        measured = _load_merfish_h5ad_slice(
            cache,
            config,
            availability["measured"],
            cell_meta.index,
            imputed=False,
        )
        if measured is not None:
            frames.append(measured)

    if availability["imputed"]:
        print(
            f"Loading {len(availability['imputed'])} imputed MERFISH genes "
            f"(use_imputed_merfish=true; first read can take several minutes)…",
            flush=True,
        )
        imputed = _load_merfish_h5ad_slice(
            cache,
            config,
            availability["imputed"],
            cell_meta.index,
            imputed=True,
        )
        if imputed is not None:
            frames.append(imputed)

    if not frames:
        return None

    if len(frames) == 1:
        combined = frames[0]
    else:
        combined = ad.concat(frames, join="outer", axis=1)
    combined.var["source"] = [
        "imputed" if g in availability["imputed"] else "measured"
        for g in combined.var["gene_symbol"]
    ]
    return combined


def aggregate_merfish_expression(
    adata: ad.AnnData,
    cell_meta: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Mean expression per (cell_type_level, brain_area) for MERFISH."""
    return aggregate_scrna_expression(adata, cell_meta, config)


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


def _join_zhuang_parcellation(
    cache: Any,
    cell: pd.DataFrame,
    dataset_id: str,
) -> pd.DataFrame:
    """Attach CCF parcellation columns to Zhuang cell metadata."""
    parc_cols = [c for c in cell.columns if c.startswith("parcellation_")]
    if parc_cols:
        return cell

    try:
        parc_cell = cache.get_metadata_dataframe(
            directory=dataset_id,
            file_name="cell_metadata_with_parcellation_annotation",
            dtype={"cell_label": str},
        )
        parc_cell = parc_cell.set_index("cell_label")
        parc_cols = [
            c for c in parc_cell.columns if c.startswith("parcellation_")
        ]
        if parc_cols:
            return cell.join(parc_cell[parc_cols], how="left")
    except Exception:
        pass

    ccf_dir = f"{dataset_id}-CCF"
    ccf_coords = cache.get_metadata_dataframe(
        directory=ccf_dir,
        file_name="ccf_coordinates",
        dtype={"cell_label": str},
    )
    ccf_coords = ccf_coords.set_index("cell_label")
    if "parcellation_index" not in ccf_coords.columns:
        warnings.warn(
            f"{ccf_dir} lacks parcellation_index; brain areas unavailable.",
            UserWarning,
            stacklevel=2,
        )
        return cell

    rename = {"x": "x_ccf", "y": "y_ccf", "z": "z_ccf"}
    ccf_coords = ccf_coords.rename(
        columns={k: v for k, v in rename.items() if k in ccf_coords.columns},
    )

    parc_ann = cache.get_metadata_dataframe(
        directory=ALLEN_CCF_DIR,
        file_name="parcellation_to_parcellation_term_membership_acronym",
    )
    parc_ann = parc_ann.set_index("parcellation_index")
    parc_ann.columns = [f"parcellation_{c}" for c in parc_ann.columns]

    cell = cell.join(
        ccf_coords[["parcellation_index"] + [c for c in rename.values() if c in ccf_coords.columns]],
        how="left",
    )
    return cell.join(parc_ann, on="parcellation_index", how="left")


def load_zhuang_cell_metadata(
    cache: Any,
    config: dict[str, Any],
    dataset_id: str,
) -> pd.DataFrame:
    """
    Zhuang MERFISH cells with taxonomy and CCF brain_area assignment.

    Filtered to config brain_areas. Indexed by cell_label.
    """
    brain_areas = config["brain_areas"]
    cell_type_col = config["cell_type_level"]

    cell = cache.get_metadata_dataframe(
        directory=dataset_id,
        file_name="cell_metadata_with_cluster_annotation",
        dtype={"cell_label": str},
    )
    cell = cell.set_index("cell_label")
    cell = _join_zhuang_parcellation(cache, cell, dataset_id)
    cell["brain_area"] = assign_merfish_brain_area(cell, brain_areas)
    cell = cell[cell["brain_area"].notna() & cell["brain_area"].isin(brain_areas)].copy()

    for ba in brain_areas:
        n = int((cell["brain_area"] == ba).sum())
        if n == 0:
            warnings.warn(
                f"No Zhuang cells in {dataset_id!r} mapped to brain_area {ba!r}",
                UserWarning,
                stacklevel=2,
            )

    if cell_type_col not in cell.columns:
        raise KeyError(
            f"cell_type_level {cell_type_col!r} not in Zhuang metadata for {dataset_id!r}",
        )

    return cell


def _zhuang_gene_panel(cache: Any, dataset_id: str) -> set[str]:
    gene_df = cache.get_metadata_dataframe(directory=dataset_id, file_name="gene")
    return set(gene_df["gene_symbol"].astype(str))


def check_zhuang_genes(
    cache: Any,
    genes: list[str],
    dataset_id: str,
) -> dict[str, list[str]]:
    """Classify genes as present or missing in the Zhuang ~1,122-gene panel."""
    panel = _zhuang_gene_panel(cache, dataset_id)
    present = [g for g in genes if g in panel]
    missing = [g for g in genes if g not in panel]
    return {"present": present, "missing": missing}


def load_zhuang_expression_subset(
    cache: Any,
    genes: list[str],
    cell_meta: pd.DataFrame,
    config: dict[str, Any],
    dataset_id: str,
) -> ad.AnnData | None:
    """Batch-load Zhuang expression for requested genes × filtered cells."""
    if not genes or cell_meta.empty:
        return None

    availability = check_zhuang_genes(cache, genes, dataset_id)
    warn_missing_genes(availability["present"], genes)

    if not availability["present"]:
        return None

    data_suffix = get_expression_suffix(config)
    path = cache.get_file_path(
        directory=dataset_id,
        file_name=f"{dataset_id}/{data_suffix}",
    )
    adata = ad.read_h5ad(path, backed="r")
    try:
        if "gene_symbol" not in adata.var.columns:
            return None

        cells_present = adata.obs_names.intersection(cell_meta.index)
        if len(cells_present) == 0:
            return None

        gene_mask = adata.var["gene_symbol"].isin(availability["present"])
        ensembl_ids = adata.var.index[gene_mask].tolist()
        if not ensembl_ids:
            return None

        subset = adata[cells_present, ensembl_ids].to_memory()
        symbol_by_ensembl = dict(
            zip(subset.var_names, subset.var["gene_symbol"].astype(str)),
        )
        subset.var["gene_symbol"] = [
            symbol_by_ensembl[e] for e in subset.var_names
        ]
        subset.var["source"] = "measured"
        return _ensure_float32_x(subset)
    finally:
        adata.file.close()


def aggregate_zhuang_replicates_mean(
    per_replicate: list[pd.DataFrame],
) -> pd.DataFrame:
    """Combine replicate-level aggregates: mean expression/detection, summed n_cells."""
    cols = [
        "cell_type", "brain_area", "gene",
        "mean_expression", "frac_expressing", "n_cells", "family",
    ]
    if not per_replicate:
        return pd.DataFrame(columns=cols)
    if len(per_replicate) == 1:
        return per_replicate[0].copy()

    combined = pd.concat(per_replicate, ignore_index=True)
    grouped = combined.groupby(
        ["cell_type", "brain_area", "gene", "family"],
        observed=True,
    )
    agg_spec: dict[str, str] = {"mean_expression": "mean"}
    if "frac_expressing" in combined.columns:
        agg_spec["frac_expressing"] = "mean"
    if "n_cells" in combined.columns:
        agg_spec["n_cells"] = "sum"
    return grouped.agg(agg_spec).reset_index()


def load_zhuang_aggregated(
    cache: Any,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Load Zhuang expression per replicate, aggregate, then mean across replicates.
    """
    datasets = get_zhuang_datasets(config)
    config["_zhuang_replicates_used"] = datasets

    per_replicate: list[pd.DataFrame] = []
    requested = list(config["_all_genes"])

    for dataset_id in tqdm(datasets, desc="Zhuang replicates"):
        availability = check_zhuang_genes(cache, requested, dataset_id)
        loadable = availability["present"]
        if not loadable:
            warnings.warn(
                f"No requested genes in Zhuang panel for {dataset_id!r}; skipping.",
                UserWarning,
                stacklevel=2,
            )
            continue

        cell_meta = load_zhuang_cell_metadata(cache, config, dataset_id)
        if cell_meta.empty:
            warnings.warn(
                f"No Zhuang cells in config brain areas for {dataset_id!r}; skipping.",
                UserWarning,
                stacklevel=2,
            )
            continue

        adata = load_zhuang_expression_subset(
            cache, loadable, cell_meta, config, dataset_id,
        )
        if adata is None:
            continue

        per_replicate.append(
            aggregate_scrna_expression(adata, cell_meta, config),
        )

    if not per_replicate:
        raise RuntimeError("No Zhuang expression data loaded for any replicate.")

    return aggregate_zhuang_replicates_mean(per_replicate)


def merfish_gene_source_map(
    cache: Any,
    genes: list[str],
    config: dict[str, Any],
) -> dict[str, str]:
    """Map gene symbol -> 'measured' or 'imputed' for Allen MERFISH."""
    availability = check_merfish_genes(cache, genes, config)
    sources: dict[str, str] = {}
    for gene in availability["measured"]:
        sources[gene] = "measured"
    for gene in availability["imputed"]:
        sources[gene] = "imputed"
    return sources


def load_allen_merfish_aggregate(
    cache: Any,
    config: dict[str, Any],
    *,
    exploration_root: Path | str | None = None,
) -> tuple[pd.DataFrame, Path | None, bool]:
    """
    Load Allen MERFISH aggregated table from a prior run parquet, or recompute.

    Returns ``(agg_long, parquet_path_or_none, reran)``.
    Emits a marked :class:`UserWarning` when re-running because no parquet was found.
    """
    dataset = config["data"].get("merfish_dataset", MERFISH_DIR)
    parquet_path = find_prior_run_parquet(
        config,
        parquet_filename="aggregated_merfish.parquet",
        dataset_slug=dataset,
        exploration_root=exploration_root,
    )

    if parquet_path is not None:
        agg = pd.read_parquet(parquet_path)
        return agg, parquet_path, False

    warnings.warn(
        "ALLEN MERFISH RE-RUN: No prior aggregated_merfish.parquet found under the "
        f"exploration folder for cell_type_level={config['cell_type_level']!r} and "
        f"dataset={dataset!r}. Re-running Allen MERFISH aggregation (may download "
        "large expression matrices). Run notebook 02 first to cache results.",
        UserWarning,
        stacklevel=2,
    )

    requested = list(config["_all_genes"])
    availability = check_merfish_genes(cache, requested, config)
    loadable = availability["measured"] + availability["imputed"]
    if not loadable:
        raise RuntimeError("No requested genes available in Allen MERFISH.")

    cell_meta = load_merfish_cell_metadata(cache, config)
    adata = load_merfish_expression_subset(cache, loadable, cell_meta, config)
    if adata is None:
        raise RuntimeError("No Allen MERFISH expression data loaded.")

    agg = aggregate_merfish_expression(adata, cell_meta, config)
    return agg, None, True


def merge_crossref_aggregates(
    allen_agg: pd.DataFrame,
    other_agg: pd.DataFrame,
    allen_gene_sources: dict[str, str],
    *,
    other_key: str = "zhuang",
) -> pd.DataFrame:
    """
    Inner-join Allen and another dataset's aggregates on cell_type × brain_area × gene.

    Adds ``allen_expression``, ``{other_key}_expression``, and ``allen_source``.
    """
    other_col = f"{other_key}_expression"
    overlap = sorted(set(allen_agg["gene"]) & set(other_agg["gene"]))
    if not overlap:
        return pd.DataFrame()

    allen = allen_agg[allen_agg["gene"].isin(overlap)].rename(
        columns={"mean_expression": "allen_expression"},
    )
    other = other_agg[other_agg["gene"].isin(overlap)].rename(
        columns={"mean_expression": other_col},
    )
    cols = ["cell_type", "brain_area", "gene", "family"]
    merged = allen[cols + ["allen_expression"]].merge(
        other[cols + [other_col]],
        on=cols,
        how="inner",
    )
    merged["allen_source"] = merged["gene"].map(allen_gene_sources).fillna("unknown")
    return merged


def _read_vizgen_gene_header(path: Path) -> list[str]:
    with open(path) as f:
        header = f.readline().strip()
    if not header:
        return []
    return [c for c in header.split(",") if c]


def _normalize_log2_cpm_plus_one(x: np.ndarray) -> np.ndarray:
    """Raw counts -> log2(CPM+1), row-wise (denominator = sum over given columns).

    Only correct when ``x`` already contains the full gene panel for each cell.
    For a gene *subset*, use :func:`_log2_cpm_plus_one_with_totals` with per-cell
    totals computed over the full panel; otherwise CPM is normalised to the subset
    and not comparable across datasets.
    """
    arr = np.asarray(x, dtype=np.float64)
    totals = arr.sum(axis=1, keepdims=True)
    return _log2_cpm_plus_one_with_totals(arr, totals)


def _log2_cpm_plus_one_with_totals(
    counts: np.ndarray,
    totals: np.ndarray,
) -> np.ndarray:
    """log2(CPM+1) for ``counts`` using externally supplied per-cell ``totals``.

    ``totals`` is the per-cell library size (sum over the full real-gene panel),
    so a gene subset is normalised against true library size rather than the
    subset sum.
    """
    arr = np.asarray(counts, dtype=np.float64)
    tot = np.asarray(totals, dtype=np.float64).reshape(-1, 1)
    tot = np.where(tot <= 0, 1.0, tot)
    cpm = arr / tot * 1e6
    return np.log2(cpm + 1.0)


def _knn_majority_labels_with_confidence(
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    x_query: np.ndarray,
    *,
    k: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Majority-vote labels among *k* nearest neighbours, with vote-fraction confidence.

    Returns ``(labels, confidence)`` where confidence is the fraction of the k
    neighbours that carried the winning label (1/k .. 1.0).
    """
    if len(x_ref) == 0:
        raise ValueError("Reference matrix is empty")

    k_eff = min(k, len(x_ref))
    tree = cKDTree(x_ref)
    _, indices = tree.query(x_query, k=k_eff)
    if k_eff == 1:
        indices = np.asarray(indices).reshape(-1, 1)

    labels: list[str] = []
    confidences: list[float] = []
    y_arr = np.asarray(y_ref)
    for row in np.asarray(indices):
        neighbours = y_arr[row]
        uniq, counts = np.unique(neighbours, return_counts=True)
        winner = counts.argmax()
        labels.append(str(uniq[winner]))
        confidences.append(float(counts[winner]) / float(k_eff))
    return np.array(labels, dtype=object), np.asarray(confidences, dtype=np.float64)


def _knn_majority_labels(
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    x_query: np.ndarray,
    *,
    k: int,
) -> np.ndarray:
    """Assign labels by majority vote among *k* nearest reference neighbours."""
    labels, _ = _knn_majority_labels_with_confidence(x_ref, y_ref, x_query, k=k)
    return labels


def _standardize_fit(x_ref: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-feature mean/std from the reference matrix (std floored to 1.0)."""
    mu = np.asarray(x_ref, dtype=np.float64).mean(axis=0)
    sigma = np.asarray(x_ref, dtype=np.float64).std(axis=0)
    sigma = np.where(sigma <= 0, 1.0, sigma)
    return mu, sigma


def _subsample_reference_cells(
    cell_meta: pd.DataFrame,
    max_cells: int,
    stratify_cols: list[str],
) -> pd.Index:
    """Stratified subsample of reference cell indices (cap at ``max_cells``)."""
    if len(cell_meta) <= max_cells:
        return cell_meta.index

    groups = cell_meta.groupby(stratify_cols, observed=True)
    n_groups = max(len(groups), 1)
    per_group = max(1, max_cells // n_groups)
    parts: list[pd.Index] = []
    for _, group in groups:
        idx = group.index
        if len(idx) <= per_group:
            parts.append(idx)
        else:
            parts.append(idx.to_series().sample(n=per_group, random_state=0).index)
    sampled = pd.Index(np.concatenate([p.to_numpy() for p in parts])).unique()
    if len(sampled) > max_cells:
        sampled = sampled.to_series().sample(n=max_cells, random_state=0).index
    return sampled


def build_allen_merfish_label_reference(
    cache: Any,
    config: dict[str, Any],
    transfer_genes: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], pd.Index]:
    """
    Subsampled Allen MERFISH expression matrix for kNN label transfer.

    Returns ``(X_ref, y_cell_type, y_brain_area, gene_order, cell_ids)``.
    """
    if not transfer_genes:
        raise ValueError("No genes available for Allen MERFISH label transfer")

    cell_meta = load_merfish_cell_metadata(cache, config)
    if cell_meta.empty:
        raise RuntimeError("No Allen MERFISH cells in config brain areas for label transfer")

    data_cfg = config.get("data", {})
    max_cells = int(data_cfg.get("vizgen_label_transfer_max_cells", 50_000))
    cell_type_col = config["cell_type_level"]
    stratify_cols = [cell_type_col, "brain_area"]

    ref_ids = _subsample_reference_cells(cell_meta, max_cells, stratify_cols)
    ref_meta = cell_meta.loc[ref_ids]

    adata = load_merfish_expression_subset(cache, transfer_genes, ref_meta, config)
    if adata is None:
        raise RuntimeError("Failed to load Allen MERFISH expression for label transfer")

    gene_order = adata.var["gene_symbol"].astype(str).tolist()
    x = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    y_type = ref_meta.loc[adata.obs_names, cell_type_col].astype(str).to_numpy()
    y_area = ref_meta.loc[adata.obs_names, "brain_area"].astype(str).to_numpy()
    return x, y_type, y_area, gene_order, adata.obs_names


def transfer_allen_merfish_labels(
    adata: ad.AnnData,
    x_ref: np.ndarray,
    y_type: np.ndarray,
    y_area: np.ndarray,
    gene_order: list[str],
    config: dict[str, Any],
) -> ad.AnnData:
    """Assign Allen ``cell_type_level`` and ``brain_area`` to Vizgen cells via kNN."""
    missing = [g for g in gene_order if g not in adata.var_names]
    if missing:
        raise ValueError(
            f"Vizgen matrix missing {len(missing)} transfer genes, e.g. {missing[:5]}",
        )

    k = int(config.get("data", {}).get("vizgen_label_transfer_k", 15))
    cell_type_col = config["cell_type_level"]

    x_query = adata[:, gene_order].X
    if hasattr(x_query, "toarray"):
        x_query = x_query.toarray()
    else:
        x_query = np.asarray(x_query)
    x_query = np.asarray(x_query, dtype=np.float64)
    x_ref = np.asarray(x_ref, dtype=np.float64)

    # Standardise per gene using reference statistics so kNN distances are not
    # dominated by high-variance genes or by cross-dataset scale differences.
    mu, sigma = _standardize_fit(x_ref)
    x_ref_z = (x_ref - mu) / sigma
    x_query_z = (x_query - mu) / sigma

    adata = adata.copy()
    type_labels, type_conf = _knn_majority_labels_with_confidence(
        x_ref_z, y_type, x_query_z, k=k,
    )
    area_labels, area_conf = _knn_majority_labels_with_confidence(
        x_ref_z, y_area, x_query_z, k=k,
    )
    adata.obs[cell_type_col] = type_labels
    adata.obs["brain_area"] = area_labels
    adata.obs["label_transfer_confidence"] = type_conf
    adata.obs["brain_area_transfer_confidence"] = area_conf
    return adata


def load_vizgen_sample(
    config: dict[str, Any],
    sample_tag: str,
    genes: list[str] | None = None,
) -> ad.AnnData:
    """Load one Vizgen replicate as AnnData (log2(CPM+1) in ``X``)."""
    data_dir = get_vizgen_data_dir(config)
    cbg_path, meta_path = vizgen_sample_file_paths(data_dir, sample_tag)

    panel = _read_vizgen_gene_header(cbg_path)
    if not panel:
        raise ValueError(f"Empty or unreadable Vizgen panel header: {cbg_path}")

    if genes:
        use_genes = [g for g in genes if g in panel]
        warn_missing_genes(use_genes, genes)
    else:
        use_genes = panel

    if not use_genes:
        raise RuntimeError(f"No requested genes found in Vizgen panel for {sample_tag}")

    all_cols = pd.read_csv(cbg_path, nrows=0).columns.tolist()
    index_col = all_cols[0]
    # Real-gene columns = everything except the index and Vizgen 'Blank*' controls.
    # Per-cell library size is computed over this full panel so that a requested
    # gene subset is normalised to true library size (CPM), not the subset sum.
    real_gene_cols = [
        c for c in all_cols
        if c != index_col and not str(c).lower().startswith("blank")
    ]
    full = pd.read_csv(cbg_path, index_col=index_col, usecols=[index_col] + real_gene_cols)
    meta = pd.read_csv(meta_path, index_col=0)

    common = full.index.intersection(meta.index)
    if len(common) == 0:
        raise RuntimeError(f"No overlapping cell IDs between expression and metadata ({sample_tag})")

    full = full.loc[common]
    meta = meta.loc[common]

    totals = full.to_numpy().sum(axis=1)
    counts = full[use_genes].to_numpy()
    x = _log2_cpm_plus_one_with_totals(counts, totals)
    adata = ad.AnnData(
        X=x.astype(np.float32),
        obs=meta,
        var=pd.DataFrame(index=use_genes),
    )
    adata.var["gene_symbol"] = use_genes
    adata.obs_names = common.astype(str)
    adata.uns["vizgen_sample"] = sample_tag
    return adata


def check_vizgen_genes(
    config: dict[str, Any],
    genes: list[str],
    sample_tag: str | None = None,
) -> dict[str, list[str]]:
    """Classify genes as present or missing in the Vizgen 483-gene panel."""
    data_dir = get_vizgen_data_dir(config)
    if sample_tag is None:
        samples = discover_vizgen_samples(config)
        if not samples:
            raise FileNotFoundError("No Vizgen samples found for gene panel check")
        sample_tag = samples[0]
    cbg_path, _ = vizgen_sample_file_paths(data_dir, sample_tag)
    panel = set(_read_vizgen_gene_header(cbg_path))
    present = [g for g in genes if g in panel]
    missing = [g for g in genes if g not in panel]
    return {"present": present, "missing": missing}


def load_vizgen_aggregated(
    cache: Any,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Load Vizgen replicates, transfer Allen labels, filter brain areas, aggregate.

    Multiple downloaded samples are averaged (mean) per cell_type × region × gene.
    """
    samples = get_vizgen_samples(config)
    config["_vizgen_samples_used"] = samples

    requested = list(config["_all_genes"])
    availability = check_vizgen_genes(config, requested, samples[0])
    loadable = availability["present"]
    if not loadable:
        raise RuntimeError("No requested genes in Vizgen panel.")

    allen_availability = check_merfish_genes(cache, requested, config)
    transfer_genes = sorted(
        set(loadable)
        & set(allen_availability["measured"] + allen_availability["imputed"]),
    )
    if not transfer_genes:
        raise RuntimeError(
            "No overlapping genes between Vizgen panel and Allen MERFISH for label transfer",
        )

    x_ref, y_type, y_area, gene_order, _ = build_allen_merfish_label_reference(
        cache, config, transfer_genes,
    )

    brain_areas = set(config["brain_areas"])
    cell_type_col = config["cell_type_level"]
    min_conf = float(
        config.get("data", {}).get("vizgen_label_transfer_min_confidence", 0.0)
    )
    per_sample: list[pd.DataFrame] = []

    for sample_tag in tqdm(samples, desc="Vizgen samples"):
        adata = load_vizgen_sample(config, sample_tag, genes=loadable)
        adata = transfer_allen_merfish_labels(
            adata, x_ref, y_type, y_area, gene_order, config,
        )
        cell_meta = adata.obs
        keep = cell_meta["brain_area"].isin(brain_areas)
        if min_conf > 0.0 and "label_transfer_confidence" in cell_meta.columns:
            keep = keep & (cell_meta["label_transfer_confidence"] >= min_conf)
        if not keep.any():
            warnings.warn(
                f"No Vizgen cells in {sample_tag!r} assigned to config brain_areas "
                f"{config['brain_areas']} after label transfer"
                + (f" with confidence ≥ {min_conf}" if min_conf > 0 else "")
                + "; skipping sample.",
                UserWarning,
                stacklevel=2,
            )
            continue

        cell_meta = cell_meta.loc[keep, [cell_type_col, "brain_area"]].copy()
        adata = adata[keep].copy()
        per_sample.append(aggregate_scrna_expression(adata, cell_meta, config))

    if not per_sample:
        raise RuntimeError(
            "No Vizgen expression aggregated for any sample; check brain_areas and "
            "label transfer settings",
        )

    return aggregate_zhuang_replicates_mean(per_sample)
