"""Functional landscape analysis: module scores, clustering, joint embedding.

Builds on the synthesis evidence table (high/medium tier only) with a
multi-level spatial expression picker so area-specific differences (e.g.
VISp vs VISpm) can use MERFISH/Vizgen/Zhuang before scRNA fallback.

Rows are (cell_type, brain_area) targets — typically supertype × region.
Experimental resonance linking is stubbed for a later milestone.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

from src.config import discover_run_parquets_by_level
from src.synthesis import DATASET_SPECS, EXPRESSED_TIERS, _dataset_slug
from src.utils import filter_cell_types_by_name

TEST_EXPERIMENTAL_RESONANCE_REL = Path("data/test_experimental_resonance.csv")

DEFAULT_EXPRESSION_PRIORITY: tuple[str, ...] = (
    "allen_merfish",
    "vizgen",
    "zhuang",
    "allen_scrna",
)

EXPERIMENTAL_RESONANCE_REQUIRED_COLUMNS: tuple[str, ...] = (
    "cell_id",
    "coarse_type",
    "brain_area_group",
    "peak_resonance_hz",
    "resonance_strength",
)
EXPERIMENTAL_RESONANCE_OPTIONAL_COLUMNS: tuple[str, ...] = (
    "brain_area",
    "subcluster_hint",
    "notes",
)
EXPERIMENTAL_RESONANCE_COLUMNS: tuple[str, ...] = (
    *EXPERIMENTAL_RESONANCE_REQUIRED_COLUMNS,
    *EXPERIMENTAL_RESONANCE_OPTIONAL_COLUMNS,
)


def functional_config(config: dict[str, Any]) -> dict[str, Any]:
    """Return the ``functional_analysis`` block with defaults filled in."""
    fa = dict(config.get("functional_analysis") or {})
    fa.setdefault("cell_type_level", config.get("cell_type_level", "supertype"))
    fa.setdefault("expression_priority", list(DEFAULT_EXPRESSION_PRIORITY))
    fa.setdefault("min_module_genes", 2)
    fa.setdefault("min_score_completeness", 0.5)
    fa.setdefault("min_targets_for_cluster", 5)
    fa.setdefault("clustering_method", "ward")
    fa.setdefault("embedding_method", "pca")
    fa.setdefault("vis_groups", {
        "VISp": ["VISp"],
        "V2M": ["VISpm", "VISam", "RSPagl"],
    })
    return fa


def default_experimental_resonance_path(project_root: Path | str) -> Path:
    """Repo-bundled schema demo CSV (``data/test_experimental_resonance.csv``)."""
    return Path(project_root).resolve() / TEST_EXPERIMENTAL_RESONANCE_REL


def resolve_experimental_resonance_path(
    config: dict[str, Any],
    project_root: Path | str,
) -> Path:
    """Configured ephys CSV path, or the repo default test fixture."""
    fa = functional_config(config)
    raw = fa.get("experimental_resonance_csv")
    if raw:
        p = Path(raw).expanduser()
        if not p.is_absolute():
            p = Path(project_root).resolve() / p
        return p
    return default_experimental_resonance_path(project_root)


def build_level_source_report(
    config: dict[str, Any],
    *,
    exploration_root: Path | str | None = None,
) -> pd.DataFrame:
    """Audit which datasets exist at each cell_type_level (from run folder names)."""
    needed = functional_config(config)["cell_type_level"]
    rows: list[dict[str, Any]] = []
    for key, spec in DATASET_SPECS.items():
        slug = _dataset_slug(config, key)
        by_level = discover_run_parquets_by_level(
            config,
            parquet_filename=spec["parquet"],
            dataset_slug=slug,
            exploration_root=exploration_root,
        )
        chosen = by_level.get(needed)
        rows.append({
            "dataset": key,
            "required_level": needed,
            "status": "ok" if chosen is not None else "MISSING",
            "run_folder": chosen.parent.name if chosen is not None else "",
            "parquet": str(chosen) if chosen is not None else "",
            "available_levels": ", ".join(sorted(by_level.keys())) or "(none)",
        })
    return pd.DataFrame(rows)


def print_level_source_report(report: pd.DataFrame) -> None:
    """Print a human-readable summary of :func:`build_level_source_report`."""
    if report.empty:
        print("No datasets in registry.")
        return
    needed = report["required_level"].iloc[0]
    print(f"Functional analysis requires cell_type_level={needed!r} run folders.")
    for _, row in report.iterrows():
        if row["status"] == "ok":
            print(f"  ✓ {row['dataset']}: {row['run_folder']}")
        else:
            print(
                f"  ✗ {row['dataset']}: MISSING at {needed!r} "
                f"(found: {row['available_levels']})"
            )
    missing = report[report["status"] == "MISSING"]
    if not missing.empty:
        warnings.warn(
            f"{len(missing)} dataset(s) missing at level {needed!r}: "
            f"{', '.join(missing['dataset'].tolist())}. "
            "Re-run the corresponding notebook(s) at that level.",
            UserWarning,
            stacklevel=2,
        )


def prepare_functional_evidence(
    config: dict[str, Any],
    *,
    exploration_root: Path | str | None = None,
    allen_gene_sources: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, dict[str, Path], pd.DataFrame]:
    """Build evidence at ``functional_analysis.cell_type_level`` from run folders."""
    from src import synthesis as syn

    level = functional_config(config)["cell_type_level"]
    report = build_level_source_report(
        config, exploration_root=exploration_root,
    )
    aggregates, sources = syn.gather_dataset_aggregates(
        config,
        exploration_root=exploration_root,
        cell_type_level=level,
    )
    if not aggregates:
        raise RuntimeError(
            f"No aggregates found at cell_type_level={level!r}. "
            "See level source report for which notebooks to re-run."
        )
    evidence = syn.build_evidence_table(
        aggregates, config, allen_gene_sources=allen_gene_sources,
    )
    return evidence, sources, report


def load_functional_modules(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Load functional module definitions from config."""
    fa = functional_config(config)
    modules = fa.get("functional_modules") or {}
    if not modules:
        raise KeyError(
            "config['functional_analysis']['functional_modules'] is empty; "
            "add resonance and neuromodulator coupling modules to query_config.yaml"
        )
    all_genes = set(config.get("_all_genes") or [])
    out: dict[str, dict[str, Any]] = {}
    for name, spec in modules.items():
        if isinstance(spec, list):
            genes = list(spec)
            meta: dict[str, Any] = {"category": "unknown"}
            raw_genes = genes
        else:
            genes = list(spec.get("genes") or [])
            meta = {k: v for k, v in spec.items() if k != "genes"}
            raw_genes = genes
        genes = [g for g in genes if g in all_genes]
        missing = [g for g in raw_genes if g not in all_genes]
        if missing:
            warnings.warn(
                f"Module {name!r}: {len(missing)} gene(s) not in panel — skipped",
                UserWarning,
                stacklevel=2,
            )
        out[name] = {**meta, "genes": genes}
    return out


def brain_area_group(region: str, config: dict[str, Any]) -> str | None:
    """Map CCF acronym to ephys grouping (VISp vs V2M)."""
    for group, regions in functional_config(config).get("vis_groups", {}).items():
        if region in regions:
            return group
    return None


def sort_targets(
    index: pd.MultiIndex,
    config: dict[str, Any],
) -> pd.MultiIndex:
    """Order ``(cell_type, brain_area)`` rows: cell type first, then brain area.

    Brain areas follow ``config['brain_areas']`` order when present.
    """
    if index.empty:
        return index
    df = index.to_frame(index=False)
    area_order = {
        str(a): i for i, a in enumerate(config.get("brain_areas") or [])
    }
    df["_area_ord"] = df["brain_area"].map(lambda a: area_order.get(str(a), 999))
    df = df.sort_values(
        ["cell_type", "_area_ord", "brain_area"],
        kind="stable",
    )
    return pd.MultiIndex.from_frame(
        df[["cell_type", "brain_area"]],
        names=["cell_type", "brain_area"],
    )


def _target_index(evidence: pd.DataFrame, config: dict[str, Any]) -> pd.MultiIndex:
    regions = set(config.get("brain_areas") or [])
    sub = evidence[evidence["brain_area"].isin(regions)].copy()
    cell_types = filter_cell_types_by_name(sub["cell_type"].unique(), config)
    sub = sub[sub["cell_type"].isin(cell_types)]
    pairs = (
        sub[["cell_type", "brain_area"]]
        .drop_duplicates()
    )
    idx = pd.MultiIndex.from_frame(pairs, names=["cell_type", "brain_area"])
    return sort_targets(idx, config)


def pick_expression_value(
    row: pd.Series,
    config: dict[str, Any],
    *,
    prefer_measured_merfish: bool = False,
) -> tuple[float, str | None, str | None]:
    """Pick one expression value using the configured dataset priority."""
    priority = functional_config(config)["expression_priority"]
    fallback: tuple[float, str | None, str | None] | None = None

    for ds in priority:
        mean_col = f"{ds}_mean"
        if mean_col not in row.index or pd.isna(row[mean_col]):
            continue
        expressed_col = f"{ds}_expressed"
        if expressed_col in row.index and not bool(row[expressed_col]):
            continue
        source_col = f"{ds}_source"
        source = row[source_col] if source_col in row.index else "measured"
        if pd.isna(source):
            source = "measured"
        value = float(row[mean_col])
        if prefer_measured_merfish and ds == "allen_merfish" and source == "imputed":
            if fallback is None:
                fallback = (value, ds, str(source))
            continue
        return value, ds, str(source)

    if fallback is not None:
        return fallback
    return np.nan, None, None


def build_target_expression_matrix(
    evidence: pd.DataFrame,
    config: dict[str, Any],
    genes: list[str] | None = None,
    *,
    confidence_tiers: frozenset[str] = EXPRESSED_TIERS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Wide expression matrix indexed by (cell_type, brain_area).

    Only high/medium confidence rows contribute values. Returns
    ``(matrix, provenance)`` where provenance holds per-target gene metadata.
    """
    genes = genes or list(config.get("_all_genes") or [])
    targets = _target_index(evidence, config)
    if targets.empty:
        raise ValueError("No (cell_type, brain_area) targets after filtering.")

    sub = evidence[
        evidence["confidence_tier"].isin(confidence_tiers)
        & evidence["gene"].isin(genes)
    ].copy()

    matrix = pd.DataFrame(index=targets, columns=genes, dtype=float)
    prov_rows: list[dict[str, Any]] = []

    for (ct, ba), gene_rows in sub.groupby(["cell_type", "brain_area"], observed=True):
        if (ct, ba) not in matrix.index:
            continue
        for _, row in gene_rows.iterrows():
            gene = row["gene"]
            prefer_meas = gene in (config.get("_gene_category") or {}) and (
                config["_gene_category"].get(gene) == "receptors"
            )
            val, ds, source = pick_expression_value(
                row, config, prefer_measured_merfish=prefer_meas,
            )
            if np.isnan(val):
                continue
            matrix.loc[(ct, ba), gene] = val
            prov_rows.append({
                "cell_type": ct,
                "brain_area": ba,
                "gene": gene,
                "expression": val,
                "dataset": ds,
                "source": source,
                "confidence_tier": row["confidence_tier"],
            })

    provenance = pd.DataFrame(prov_rows)
    matrix = matrix.loc[sort_targets(matrix.index, config)]
    return matrix, provenance


def zscore_across_targets(matrix: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene column across (cell_type, brain_area) targets."""
    out = matrix.astype(float).copy()
    for col in out.columns:
        s = out[col]
        valid = s.dropna()
        if len(valid) < 2:
            out[col] = np.nan
            continue
        mu, sigma = valid.mean(), valid.std(ddof=0)
        if sigma == 0 or np.isnan(sigma):
            out[col] = 0.0
        else:
            out[col] = (s - mu) / sigma
    return out


def compute_module_scores(
    matrix: pd.DataFrame,
    modules: dict[str, dict[str, Any]],
    *,
    min_genes: int | None = None,
) -> pd.DataFrame:
    """Mean z-scored expression per module per target."""
    z = zscore_across_targets(matrix)
    scores: dict[str, pd.Series] = {}
    for name, spec in modules.items():
        genes = spec["genes"]
        present = [g for g in genes if g in z.columns]
        if not present:
            scores[name] = pd.Series(np.nan, index=matrix.index)
            continue
        sub = z[present]
        need = min_genes if min_genes is not None else min(2, len(present))
        valid = sub.notna().sum(axis=1) >= need
        s = sub.mean(axis=1, skipna=True)
        s[~valid] = np.nan
        scores[name] = s
    return pd.DataFrame(scores)


def _valid_targets(scores: pd.DataFrame, min_frac: float = 0.5) -> pd.Index:
    """Targets with at least ``min_frac`` of module scores defined."""
    ok = scores.notna().mean(axis=1) >= min_frac
    return scores.index[ok]


def diagnose_module_scores(
    matrix: pd.DataFrame,
    scores: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Summarise why clustering may skip targets (tier sparsity, module gaps)."""
    fa = functional_config(config)
    min_frac = float(fa.get("min_score_completeness", 0.5))
    n_modules = len(scores.columns)
    need_modules = int(np.ceil(min_frac * n_modules)) if n_modules else 0
    genes_per_target = matrix.notna().sum(axis=1)
    modules_per_target = scores.notna().sum(axis=1)
    valid = _valid_targets(scores, min_frac=min_frac)
    tier_sub = matrix.index  # targets already tier-filtered at build time
    return {
        "cell_type_level": fa["cell_type_level"],
        "n_targets": len(scores),
        "n_modules": n_modules,
        "min_score_completeness": min_frac,
        "min_modules_required": need_modules,
        "n_targets_clusterable": len(valid),
        "min_targets_for_cluster": fa["min_targets_for_cluster"],
        "genes_per_target_median": float(genes_per_target.median()) if len(genes_per_target) else 0,
        "genes_per_target_max": int(genes_per_target.max()) if len(genes_per_target) else 0,
        "modules_per_target_median": float(modules_per_target.median()) if len(modules_per_target) else 0,
        "modules_per_target_max": int(modules_per_target.max()) if len(modules_per_target) else 0,
        "n_provenance_rows": None,  # filled by caller if desired
    }


def cluster_targets(
    scores: pd.DataFrame,
    config: dict[str, Any],
    *,
    min_frac: float | None = None,
) -> pd.Series:
    """Hierarchical clustering on module scores; returns cluster labels."""
    fa = functional_config(config)
    if min_frac is None:
        min_frac = float(fa.get("min_score_completeness", 0.5))
    valid = _valid_targets(scores, min_frac=min_frac)
    if len(valid) < fa["min_targets_for_cluster"]:
        warnings.warn(
            f"Only {len(valid)} targets pass module-score completeness filter; "
            f"need {fa['min_targets_for_cluster']} for clustering.",
            UserWarning,
            stacklevel=2,
        )
        return pd.Series(dtype=int)

    sub = scores.loc[valid].fillna(0.0)
    method = fa["clustering_method"]
    dist = pdist(sub.to_numpy(), metric="euclidean")
    link = hierarchy.linkage(dist, method=method)
    n_clusters = max(2, min(8, len(valid) // 3))
    labels = hierarchy.fcluster(link, t=n_clusters, criterion="maxclust")
    return pd.Series(labels, index=valid, name="cluster")


def joint_embedding(
    scores: pd.DataFrame,
    config: dict[str, Any],
    *,
    min_frac: float | None = None,
    n_components: int = 2,
) -> pd.DataFrame:
    """PCA (default) or UMAP on module scores for joint resonance × neuromod view."""
    fa = functional_config(config)
    if min_frac is None:
        min_frac = float(fa.get("min_score_completeness", 0.5))
    valid = _valid_targets(scores, min_frac=min_frac)
    if len(valid) < 3:
        return pd.DataFrame(index=scores.index)

    sub = scores.loc[valid].fillna(0.0)
    method = str(fa.get("embedding_method", "pca")).lower()

    if method == "umap":
        try:
            import umap  # type: ignore[import-untyped]
        except ImportError:
            warnings.warn(
                "umap-learn not installed; falling back to PCA for joint embedding.",
                UserWarning,
                stacklevel=2,
            )
            method = "pca"

    if method == "umap":
        reducer = umap.UMAP(n_components=n_components, random_state=0)
        coords = reducer.fit_transform(sub.to_numpy())
        cols = [f"UMAP{i + 1}" for i in range(n_components)]
    else:
        x = sub.to_numpy()
        x = x - x.mean(axis=0)
        u, s, _vt = np.linalg.svd(x, full_matrices=False)
        coords = u[:, :n_components] * s[:n_components]
        cols = [f"PC{i + 1}" for i in range(n_components)]
        var = (s ** 2) / max(len(valid) - 1, 1)
        total = var.sum() or 1.0
        cols = [
            f"{c} ({100 * var[i] / total:.0f}% var)" if i < len(var) else c
            for i, c in enumerate(cols)
        ]

    emb = pd.DataFrame(coords, index=valid, columns=cols[:n_components])
    return emb.reindex(scores.index)


def target_labels(index: pd.MultiIndex) -> list[str]:
    return [f"{ct} × {ba}" for ct, ba in index]


def generate_test_experimental_resonance(
    path: Path | str | None = None,
    *,
    cells_per_type: int = 20,
    seed: int = 42,
) -> pd.DataFrame:
    """Synthetic **single-cell** resonance table (n cells per coarse_type × area).

    L5 ET / VISp is bimodal (two subclusters) to demo within-type heterogeneity
    for a future expression↔ephys clustering milestone.
    """
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []

    def _add_cells(
        prefix: str,
        coarse: str,
        bag: str,
        n: int,
        peak_mu: float,
        peak_sigma: float,
        strength_mu: float,
        strength_sigma: float,
        *,
        subcluster_hint: str | None = None,
        brain_area: str | None = None,
        missing_peak: int = 0,
        notes: str = "",
    ) -> None:
        for i in range(n):
            cell_id = f"{prefix}_{i + 1:02d}"
            peak = (
                float(rng.normal(peak_mu, peak_sigma))
                if i >= missing_peak
                else np.nan
            )
            strength = (
                float(rng.normal(strength_mu, strength_sigma))
                if i >= missing_peak
                else np.nan
            )
            rows.append({
                "cell_id": cell_id,
                "coarse_type": coarse,
                "brain_area_group": bag,
                "brain_area": brain_area or "",
                "peak_resonance_hz": peak,
                "resonance_strength": max(0.05, strength) if pd.notna(strength) else np.nan,
                "subcluster_hint": subcluster_hint or "",
                "notes": notes if i == 0 and notes else "",
            })

    n = cells_per_type
    half = n // 2

    # L5 IT — unimodal, theta-ish M-resonance
    _add_cells("L5IT_VISp", "L5 IT", "VISp", n, 2.4, 0.35, 0.48, 0.07, brain_area="VISp")
    _add_cells("L5IT_V2M", "L5 IT", "V2M", n, 1.9, 0.40, 0.36, 0.08,
               missing_peak=2, notes="2 cells missing peak (QC fail)")

    # L5 ET — bimodal in VISp (within-type subclusters); unimodal in V2M
    _add_cells(
        "L5ET_VISp_A", "L5 ET", "VISp", half, 3.9, 0.30, 0.52, 0.06,
        subcluster_hint="ET_low", brain_area="VISp",
        notes="L5 ET subcluster A — lower peak (hypothesized functional subgroup)",
    )
    _add_cells(
        "L5ET_VISp_B", "L5 ET", "VISp", n - half, 6.3, 0.35, 0.64, 0.07,
        subcluster_hint="ET_high", brain_area="VISp",
        notes="L5 ET subcluster B — higher peak",
    )
    _add_cells("L5ET_V2M", "L5 ET", "V2M", n, 4.5, 0.45, 0.55, 0.08, brain_area="VISpm")

    # L2/3 IT reference
    _add_cells("L23IT_VISp", "L2/3 IT", "VISp", n, 3.1, 0.30, 0.42, 0.06, brain_area="VISp")

    # Edge case: unmapped area group
    rows.append({
        "cell_id": "L5ET_UNKNOWN_01",
        "coarse_type": "L5 ET",
        "brain_area_group": "UNKNOWN",
        "brain_area": "",
        "peak_resonance_hz": 5.8,
        "resonance_strength": 0.58,
        "subcluster_hint": "",
        "notes": "edge case: unmapped brain_area_group",
    })

    df = pd.DataFrame(rows, columns=list(EXPERIMENTAL_RESONANCE_COLUMNS))
    if path is not None:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False)
    return df


def load_experimental_resonance(path: Path | str | None) -> pd.DataFrame:
    """Load single-cell experimental CSV; return empty typed frame if missing."""
    cols = list(EXPERIMENTAL_RESONANCE_COLUMNS)
    if path is None or not Path(path).is_file():
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path)
    missing = set(EXPERIMENTAL_RESONANCE_REQUIRED_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(
            f"Experimental resonance CSV missing required columns: {sorted(missing)}"
        )
    for opt in EXPERIMENTAL_RESONANCE_OPTIONAL_COLUMNS:
        if opt not in df.columns:
            df[opt] = ""
    if df["cell_id"].duplicated().any():
        dupes = df.loc[df["cell_id"].duplicated(), "cell_id"].tolist()
        raise ValueError(f"Duplicate cell_id values: {dupes[:5]}")
    return df[cols]


def summarize_experimental_match(
    experimental: pd.DataFrame,
    coarse_type: str | None,
    brain_area_group: str | None,
) -> dict[str, Any]:
    """Aggregate single-cell ephys stats for a coarse_type × area group."""
    empty: dict[str, Any] = {
        "ephys_n_cells": 0,
        "ephys_n_with_peak": 0,
        "ephys_peak_resonance_hz_mean": np.nan,
        "ephys_peak_resonance_hz_std": np.nan,
        "ephys_resonance_strength_mean": np.nan,
        "ephys_subcluster_hints": "",
        "ephys_match": "none",
    }
    if experimental.empty or not coarse_type or not brain_area_group:
        return empty

    match = experimental[
        (experimental["coarse_type"] == coarse_type)
        & (experimental["brain_area_group"] == brain_area_group)
    ]
    if match.empty:
        return empty

    valid = match.dropna(subset=["peak_resonance_hz"])
    hints = sorted(h for h in match["subcluster_hint"].dropna().unique() if str(h).strip())
    out = {
        "ephys_n_cells": len(match),
        "ephys_n_with_peak": len(valid),
        "ephys_peak_resonance_hz_mean": float(valid["peak_resonance_hz"].mean()) if len(valid) else np.nan,
        "ephys_peak_resonance_hz_std": float(valid["peak_resonance_hz"].std(ddof=1)) if len(valid) > 1 else np.nan,
        "ephys_resonance_strength_mean": float(valid["resonance_strength"].mean()) if len(valid) else np.nan,
        "ephys_subcluster_hints": ", ".join(hints),
        "ephys_match": "cell_level_summary",
    }
    return out


def experimental_within_type_summary(
    experimental: pd.DataFrame,
) -> pd.DataFrame:
    """Distribution of single-cell resonance per coarse_type × brain_area_group.

    Highlights within-type spread (e.g. bimodal L5 ET) for future clustering
    against expression module scores at cell or pseudobulk resolution.
    """
    if experimental.empty:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    group_cols = ["coarse_type", "brain_area_group"]
    for keys, grp in experimental.groupby(group_cols, observed=True):
        coarse, bag = keys
        valid = grp.dropna(subset=["peak_resonance_hz"])
        peaks = valid["peak_resonance_hz"]
        hints = grp["subcluster_hint"].dropna().astype(str).str.strip()
        hints = sorted(h for h in hints.unique() if h)
        rows.append({
            "coarse_type": coarse,
            "brain_area_group": bag,
            "n_cells": len(grp),
            "n_with_peak": len(valid),
            "peak_hz_mean": peaks.mean() if len(valid) else np.nan,
            "peak_hz_std": peaks.std(ddof=1) if len(valid) > 1 else np.nan,
            "peak_hz_min": peaks.min() if len(valid) else np.nan,
            "peak_hz_max": peaks.max() if len(valid) else np.nan,
            "strength_mean": valid["resonance_strength"].mean() if len(valid) else np.nan,
            "subcluster_hints": ", ".join(hints),
            "suggestive_bimodal": len(hints) >= 2,
        })
    return pd.DataFrame(rows).sort_values(group_cols)


def infer_coarse_type(cell_type: str) -> str | None:
    """Map Allen cell type name to coarse ephys bucket (best-effort substring)."""
    name = cell_type.upper()
    if "L2/3" in cell_type or "L23" in name:
        return "L2/3 IT"
    if "L5 ET" in cell_type or "L5ET" in name.replace(" ", ""):
        return "L5 ET"
    if "L5 IT" in cell_type or "L5/6 IT" in cell_type:
        return "L5 IT"
    if "L4/5 IT" in cell_type:
        return "L5 IT"  # approximate for ephys comparison
    return None


def compute_vis_group_contrast(
    module_scores: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Per cell type: V2M mean − VISp mean for each module score."""
    module_cols = list(module_scores.columns)
    long_rows: list[dict[str, Any]] = []
    for (ct, ba), score_row in module_scores.iterrows():
        grp = brain_area_group(str(ba), config)
        if grp is None:
            continue
        rec: dict[str, Any] = {
            "cell_type": ct,
            "brain_area": ba,
            "vis_group": grp,
            "coarse_type": infer_coarse_type(str(ct)),
        }
        for col in module_cols:
            rec[col] = score_row[col]
        long_rows.append(rec)

    if not long_rows:
        return pd.DataFrame()

    long = pd.DataFrame(long_rows)
    contrast_rows: list[dict[str, Any]] = []
    for ct, grp_df in long.groupby("cell_type", observed=True):
        visp = grp_df[grp_df["vis_group"] == "VISp"]
        v2m = grp_df[grp_df["vis_group"] == "V2M"]
        if visp.empty or v2m.empty:
            continue
        visp_mean = visp[module_cols].mean(numeric_only=True)
        v2m_mean = v2m[module_cols].mean(numeric_only=True)
        delta = v2m_mean - visp_mean
        rec = {
            "cell_type": ct,
            "coarse_type": infer_coarse_type(str(ct)),
            "n_VISp_regions": len(visp),
            "n_V2M_regions": len(v2m),
        }
        for col in module_cols:
            rec[f"VISp_{col}"] = visp_mean[col]
            rec[f"V2M_{col}"] = v2m_mean[col]
            rec[f"delta_{col}"] = delta[col]
        contrast_rows.append(rec)

    out = pd.DataFrame(contrast_rows)
    if not out.empty:
        out = out.sort_values(["coarse_type", "cell_type"], kind="stable")
    return out


def link_targets_to_experimental(
    module_scores: pd.DataFrame,
    experimental: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Placeholder join: atlas targets ↔ single-cell ephys cohort summaries.

    Each expression target (supertype × region) is linked to **all** ephys cells
    sharing its coarse_type and brain_area_group (VISp vs V2M). Per-target rows
    carry module scores plus aggregated ephys stats (mean/std peak Hz, n cells).

    Future milestone: cluster ephys cells within coarse types (e.g. L5 ET) and
    correlate subclusters with expression module scores at matched resolution.
    """
    rows: list[dict[str, Any]] = []
    for (ct, ba), score_row in module_scores.iterrows():
        bag = brain_area_group(str(ba), config)
        coarse = infer_coarse_type(str(ct))
        rec: dict[str, Any] = {
            "cell_type": ct,
            "brain_area": ba,
            "brain_area_group": bag,
            "coarse_type": coarse,
            "target_label": f"{ct} × {ba}",
        }
        for col in module_scores.columns:
            rec[f"score_{col}"] = score_row[col]

        rec.update(summarize_experimental_match(experimental, coarse, bag))
        rows.append(rec)

    return pd.DataFrame(rows)


def plot_module_score_heatmap(
    scores: pd.DataFrame,
    config: dict[str, Any],
    *,
    output_path: Path | str | None = None,
    title: str = "Functional module scores",
) -> plt.Figure:
    """Heatmap of module scores (targets × modules)."""
    fa_cfg = config.get("functional_analysis", {}) or {}
    cmap = config.get("output", {}).get("heatmap_cmap", "RdBu_r")
    figsize = tuple(config.get("output", {}).get("figsize_heatmap", [14, 8]))
    fig, ax = plt.subplots(figsize=figsize)

    plot_df = scores.loc[sort_targets(scores.index, config)].copy()
    plot_df.index = target_labels(plot_df.index)
    sns.heatmap(
        plot_df,
        cmap=cmap,
        center=0,
        ax=ax,
        cbar_kws={"label": "module score (z-mean)"},
        linewidths=0.2,
    )
    ax.set_title(title)
    ax.set_xlabel("Module")
    ax.set_ylabel("Target (cell_type × brain_area)")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=config.get("output", {}).get("dpi", 150), bbox_inches="tight")
    return fig


def plot_joint_embedding(
    embedding: pd.DataFrame,
    module_scores: pd.DataFrame,
    config: dict[str, Any],
    *,
    color_module: str | None = None,
    output_path: Path | str | None = None,
    title: str = "Joint functional embedding",
) -> plt.Figure:
    """Scatter of PCA/UMAP coords; colour by module score or brain area group."""
    if embedding.shape[1] < 2:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "Not enough targets for 2D embedding", ha="center")
        return fig

    xcol, ycol = embedding.columns[:2]
    fig, ax = plt.subplots(figsize=(8, 7))

    plot_idx = embedding.dropna(how="all").index
    x = embedding.loc[plot_idx, xcol]
    y = embedding.loc[plot_idx, ycol]

    if color_module and color_module in module_scores.columns:
        c = module_scores.loc[plot_idx, color_module]
        sc = ax.scatter(x, y, c=c, cmap="viridis", s=60, edgecolors="k", linewidths=0.3)
        plt.colorbar(sc, ax=ax, label=color_module)
    else:
        groups = [brain_area_group(ba, config) or ba for _, ba in plot_idx]
        palette = dict(zip(sorted(set(groups)), sns.color_palette("tab10", n_colors=len(set(groups)))))
        for g in sorted(set(groups)):
            mask = [gr == g for gr in groups]
            ax.scatter(x[mask], y[mask], label=g, s=60, edgecolors="k", linewidths=0.3)
        ax.legend(title="Area group", loc="best", fontsize=8)

    for idx, (xi, yi) in zip(plot_idx, zip(x, y)):
        ct, ba = idx
        short = ct.split(" ", 1)[-1][:28] if " " in ct else ct[:28]
        ax.annotate(f"{short}\n{ba}", (xi, yi), fontsize=5, alpha=0.75, ha="center")

    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(title)
    ax.axhline(0, color="0.85", lw=0.8)
    ax.axvline(0, color="0.85", lw=0.8)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=config.get("output", {}).get("dpi", 150), bbox_inches="tight")
    return fig


def plot_cluster_dendrogram(
    scores: pd.DataFrame,
    config: dict[str, Any],
    *,
    output_path: Path | str | None = None,
    title: str = "Hierarchical clustering on module scores",
) -> plt.Figure:
    """Dendrogram from module-score Euclidean distance."""
    min_frac = float(functional_config(config).get("min_score_completeness", 0.5))
    valid = sort_targets(_valid_targets(scores, min_frac=min_frac), config)
    sub = scores.loc[valid].fillna(0.0)
    fig, ax = plt.subplots(figsize=(12, 5))
    if len(valid) < 3:
        ax.text(0.5, 0.5, "Too few targets for dendrogram", ha="center")
        return fig

    dist = pdist(sub.to_numpy(), metric="euclidean")
    link = hierarchy.linkage(dist, method=functional_config(config)["clustering_method"])
    labels = [f"{ct} | {ba}" for ct, ba in valid]
    hierarchy.dendrogram(link, labels=labels, ax=ax, leaf_rotation=90, leaf_font_size=7)
    ax.set_title(title)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=config.get("output", {}).get("dpi", 150), bbox_inches="tight")
    return fig


def export_functional_landscape(
    matrix: pd.DataFrame,
    scores: pd.DataFrame,
    embedding: pd.DataFrame,
    clusters: pd.Series,
    provenance: pd.DataFrame,
    link_table: pd.DataFrame,
    output_dir: Path | str,
    *,
    vis_contrast: pd.DataFrame | None = None,
    level_report: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """Write functional analysis outputs to ``output_dir / functional/``."""
    out = Path(output_dir) / "functional"
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    mat_path = out / "expression_matrix.parquet"
    matrix = matrix.loc[sort_targets(matrix.index, config)]
    matrix.reset_index().to_parquet(mat_path, index=False)
    paths["expression_matrix"] = mat_path

    score_path = out / "module_scores.parquet"
    scores = scores.loc[sort_targets(scores.index, config)]
    scores.reset_index().to_parquet(score_path, index=False)
    paths["module_scores"] = score_path

    if not embedding.empty:
        emb_path = out / "joint_embedding.parquet"
        embedding.reset_index().to_parquet(emb_path, index=False)
        paths["joint_embedding"] = emb_path

    if not clusters.empty:
        cl_path = out / "clusters.csv"
        clusters.reset_index(name="cluster").to_csv(cl_path, index=False)
        paths["clusters"] = cl_path

    prov_path = out / "expression_provenance.parquet"
    provenance.to_parquet(prov_path, index=False)
    paths["provenance"] = prov_path

    link_path = out / "experimental_link_placeholder.csv"
    link_table.to_csv(link_path, index=False)
    paths["experimental_link"] = link_path

    if vis_contrast is not None and not vis_contrast.empty:
        vc_path = out / "vis_group_contrast.csv"
        vis_contrast.to_csv(vc_path, index=False)
        paths["vis_group_contrast"] = vc_path

    if level_report is not None and not level_report.empty:
        lr_path = out / "level_source_report.csv"
        level_report.to_csv(lr_path, index=False)
        paths["level_source_report"] = lr_path

    return paths


def export_experimental_sidecar(
    experimental: pd.DataFrame,
    output_dir: Path | str,
) -> dict[str, Path]:
    """Copy ephys cell table + within-type distribution summary to functional/."""
    out = Path(output_dir) / "functional"
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    if not experimental.empty:
        cell_path = out / "experimental_cells.csv"
        experimental.to_csv(cell_path, index=False)
        paths["experimental_cells"] = cell_path
        summary = experimental_within_type_summary(experimental)
        if not summary.empty:
            sum_path = out / "experimental_within_type_summary.csv"
            summary.to_csv(sum_path, index=False)
            paths["experimental_within_type_summary"] = sum_path
    return paths
