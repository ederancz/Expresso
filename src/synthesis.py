"""Cross-dataset synthesis: per-gene evidence tables, concordance, confidence tiers.

Level-agnostic. Consumes the per-dataset aggregated long tables produced by
notebooks 01-04 (``aggregated_scrna/merfish/vizgen/zhuang.parquet``) and builds:

- a tidy **evidence table** keyed by (cell_type, brain_area, gene) with per-dataset
  mean expression, detection rate, detection flag and provenance;
- **cross-dataset concordance** counts (independent measured datasets only);
- a per-row **confidence tier** (high / medium / low);
- focused dot plots and a textual statement scaffold for any (cell_type, region).

Independence model
------------------
Allen MERFISH *imputed* values are predicted from WMB-10x scRNA, so they are not
independent of Allen scRNA. They are tracked as *supporting* evidence only.
Independent measured datasets: Allen scRNA, Allen MERFISH (measured genes),
Vizgen (measured), Zhuang (measured).
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.config import find_prior_run_parquet, resolve_output_dir

# Dataset registry. ``slug`` matches the ``dataset`` passed to start_run in each
# notebook; ``parquet`` is the per-run aggregate filename.
DATASET_SPECS: dict[str, dict[str, Any]] = {
    "allen_scrna": {
        "parquet": "aggregated_scrna.parquet",
        "slug": "WMB-10Xv3",
        "independent": True,
        "region_resolved": False,
        "has_imputed": False,
        "label": "Allen scRNA",
    },
    "allen_merfish": {
        "parquet": "aggregated_merfish.parquet",
        "slug": None,  # filled from config['data']['merfish_dataset']
        "independent": True,
        "region_resolved": True,
        "has_imputed": True,
        "label": "Allen MERFISH",
    },
    "vizgen": {
        "parquet": "aggregated_vizgen.parquet",
        "slug": "Vizgen-MERFISH",
        "independent": True,
        "region_resolved": True,
        "has_imputed": False,
        "label": "Vizgen MERFISH",
    },
    "zhuang": {
        "parquet": "aggregated_zhuang.parquet",
        "slug": "Zhuang-ABCA",
        "independent": True,
        "region_resolved": True,
        "has_imputed": False,
        "label": "Zhuang MERFISH",
    },
}

DEFAULT_MIN_FRAC = 0.25
DEFAULT_MIN_MEAN = 1.0


def _dataset_slug(config: dict[str, Any], dataset_key: str) -> str:
    spec = DATASET_SPECS[dataset_key]
    if spec["slug"] is not None:
        return spec["slug"]
    if dataset_key == "allen_merfish":
        return config["data"].get("merfish_dataset", "MERFISH-C57BL6J-638850")
    raise KeyError(dataset_key)


def discover_dataset_parquets(
    config: dict[str, Any],
    *,
    exploration_root: Path | str | None = None,
) -> dict[str, Path]:
    """Locate the newest matching per-dataset aggregate parquet for each dataset."""
    found: dict[str, Path] = {}
    for key, spec in DATASET_SPECS.items():
        path = find_prior_run_parquet(
            config,
            parquet_filename=spec["parquet"],
            dataset_slug=_dataset_slug(config, key),
            exploration_root=exploration_root,
        )
        if path is not None:
            found[key] = path
    return found


def gather_dataset_aggregates(
    config: dict[str, Any],
    *,
    exploration_root: Path | str | None = None,
) -> tuple[dict[str, pd.DataFrame], dict[str, Path]]:
    """Read available per-dataset aggregate parquets.

    Returns ``(aggregates, sources)`` mapping dataset_key -> DataFrame / parquet path.
    Datasets without a discoverable parquet are skipped with a warning (re-run the
    matching notebook to produce them).
    """
    parquets = discover_dataset_parquets(config, exploration_root=exploration_root)
    aggregates: dict[str, pd.DataFrame] = {}
    for key in DATASET_SPECS:
        path = parquets.get(key)
        if path is None:
            warnings.warn(
                f"No aggregate parquet found for {key!r} "
                f"(expected {DATASET_SPECS[key]['parquet']} under a "
                f"{config['cell_type_level']}/{_dataset_slug(config, key)} run). "
                "Run the matching notebook first.",
                UserWarning,
                stacklevel=2,
            )
            continue
        aggregates[key] = pd.read_parquet(path)
    return aggregates, parquets


def _ensure_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Guarantee mean_expression / frac_expressing / n_cells exist."""
    out = df.copy()
    if "mean_expression" not in out.columns:
        raise KeyError("aggregate is missing 'mean_expression'")
    if "frac_expressing" not in out.columns:
        out["frac_expressing"] = np.nan
    if "n_cells" not in out.columns:
        out["n_cells"] = np.nan
    return out


def _detection_flag(
    mean: pd.Series,
    frac: pd.Series,
    *,
    min_frac: float,
    min_mean: float,
) -> pd.Series:
    """Expressed if mean ≥ min_mean AND (frac ≥ min_frac when frac available)."""
    by_mean = mean >= min_mean
    by_frac = frac.isna() | (frac >= min_frac)
    return by_mean & by_frac


def build_evidence_table(
    aggregates: dict[str, pd.DataFrame],
    config: dict[str, Any],
    *,
    allen_gene_sources: dict[str, str] | None = None,
    min_frac: float | None = None,
    min_mean: float | None = None,
) -> pd.DataFrame:
    """Build the per-(cell_type, brain_area, gene) cross-dataset evidence table.

    Region-resolved datasets (Allen/Vizgen/Zhuang MERFISH) are merged on
    (cell_type, brain_area, gene). Allen scRNA (not region-resolved) is merged on
    (cell_type, gene), broadcast across regions and flagged ``allen_scrna_region_resolved``.
    """
    synth_cfg = config.get("synthesis", {}) or {}
    min_frac = float(synth_cfg.get("min_frac", DEFAULT_MIN_FRAC)) if min_frac is None else min_frac
    min_mean = float(synth_cfg.get("min_mean", DEFAULT_MIN_MEAN)) if min_mean is None else min_mean
    allen_gene_sources = allen_gene_sources or {}
    genes_flat = config.get("_genes_flat", {})
    gene_category = config.get("_gene_category", {}) or {}

    region_keys = [k for k in aggregates if DATASET_SPECS[k]["region_resolved"]]
    if not region_keys:
        raise RuntimeError(
            "No region-resolved datasets available; need at least one of "
            "allen_merfish / vizgen / zhuang to anchor (cell_type, region, gene)."
        )

    # Master key set from region-resolved datasets.
    master_parts = []
    for key in region_keys:
        df = _ensure_columns(aggregates[key])
        master_parts.append(df[["cell_type", "brain_area", "gene"]])
    master = (
        pd.concat(master_parts, ignore_index=True)
        .drop_duplicates()
        .reset_index(drop=True)
    )
    master["family"] = master["gene"].map(genes_flat).fillna("unknown")
    if gene_category:
        master["category"] = master["gene"].map(gene_category).fillna("unknown")

    for key in DATASET_SPECS:
        if key not in aggregates:
            continue
        spec = DATASET_SPECS[key]
        df = _ensure_columns(aggregates[key])

        master[f"{key}_region_resolved"] = spec["region_resolved"]
        if spec["region_resolved"]:
            sub = df[["cell_type", "brain_area", "gene", "mean_expression",
                      "frac_expressing", "n_cells"]].copy()
            merged = master.merge(
                sub, on=["cell_type", "brain_area", "gene"], how="left",
            )
        else:
            # Collapse to (cell_type, gene) across pooled regions, broadcast.
            grp = (
                df.groupby(["cell_type", "gene"], observed=True)
                .agg(mean_expression=("mean_expression", "mean"),
                     frac_expressing=("frac_expressing", "mean"),
                     n_cells=("n_cells", "sum"))
                .reset_index()
            )
            merged = master.merge(grp, on=["cell_type", "gene"], how="left")

        master[f"{key}_mean"] = merged["mean_expression"].to_numpy()
        master[f"{key}_frac"] = merged["frac_expressing"].to_numpy()
        master[f"{key}_n_cells"] = merged["n_cells"].to_numpy()

        # Provenance per gene (measured vs imputed) — only Allen MERFISH varies.
        if spec["has_imputed"]:
            master[f"{key}_source"] = master["gene"].map(
                lambda g: allen_gene_sources.get(g, "measured")
            )
        else:
            present = merged["mean_expression"].notna()
            master[f"{key}_source"] = [
                "measured" if p else None for p in present
            ]

        master[f"{key}_expressed"] = _detection_flag(
            master[f"{key}_mean"],
            master[f"{key}_frac"],
            min_frac=min_frac,
            min_mean=min_mean,
        ).fillna(False)

    # Concordance: independent measured detections vs supporting imputed.
    indep_measured = pd.Series(0, index=master.index, dtype=int)
    n_present = pd.Series(0, index=master.index, dtype=int)
    supporting_imputed = pd.Series(False, index=master.index)

    for key in DATASET_SPECS:
        if key not in aggregates:
            continue
        spec = DATASET_SPECS[key]
        present = master[f"{key}_mean"].notna()
        n_present = n_present + present.astype(int)
        expressed = master[f"{key}_expressed"].fillna(False)
        is_measured = master[f"{key}_source"] == "measured"
        if spec["independent"]:
            indep_measured = indep_measured + (expressed & is_measured).astype(int)
        if spec["has_imputed"]:
            supporting_imputed = supporting_imputed | (
                expressed & (master[f"{key}_source"] == "imputed")
            )

    master["n_datasets_present"] = n_present
    master["n_independent_measured_detections"] = indep_measured
    master["supporting_imputed_detection"] = supporting_imputed

    def _tier(row: pd.Series) -> str:
        if row["n_independent_measured_detections"] >= 2:
            return "high"
        if row["n_independent_measured_detections"] == 1:
            return "medium"
        if row["supporting_imputed_detection"]:
            return "low"
        return "none"

    master["confidence_tier"] = master.apply(_tier, axis=1)
    return master


def target_evidence(
    evidence: pd.DataFrame,
    *,
    cell_type: str,
    region: str | None = None,
) -> pd.DataFrame:
    """Slice the evidence table to a cell type (and optional region), sorted by tier."""
    sub = evidence[evidence["cell_type"] == cell_type]
    if region is not None:
        sub = sub[sub["brain_area"] == region]
    tier_order = {"high": 0, "medium": 1, "low": 2, "none": 3}
    return sub.assign(_o=sub["confidence_tier"].map(tier_order)).sort_values(
        ["_o", "family", "gene"]
    ).drop(columns="_o").reset_index(drop=True)


def summarize_target(
    evidence: pd.DataFrame,
    *,
    cell_type: str,
    region: str | None = None,
) -> dict[str, Any]:
    """Return tier->gene lists and per-family expressed genes for a target."""
    sub = target_evidence(evidence, cell_type=cell_type, region=region)
    tiers = {
        t: sorted(sub.loc[sub["confidence_tier"] == t, "gene"].tolist())
        for t in ("high", "medium", "low")
    }
    expressed = sub[sub["confidence_tier"].isin(["high", "medium"])]
    families = {
        fam: sorted(g["gene"].tolist())
        for fam, g in expressed.groupby("family", observed=True)
    }
    families_by_category: dict[str, dict[str, list[str]]] = {}
    if "category" in expressed.columns:
        for cat, cdf in expressed.groupby("category", observed=True):
            families_by_category[str(cat)] = {
                fam: sorted(g["gene"].tolist())
                for fam, g in cdf.groupby("family", observed=True)
            }
    return {
        "cell_type": cell_type,
        "region": region,
        "tiers": tiers,
        "families_expressed": families,
        "families_by_category": families_by_category,
        "n_genes_considered": int(len(sub)),
    }


def statement_scaffold(summary: dict[str, Any]) -> str:
    """Human-readable scaffold for the cross-validated expression claim."""
    ct = summary["cell_type"]
    region = summary["region"] or "(region-pooled)"
    high = summary["tiers"]["high"]
    medium = summary["tiers"]["medium"]
    by_cat = summary.get("families_by_category") or {}
    lines = [
        f"In {region}, {ct} neurons express (cross-validated):",
        f"  High confidence (≥2 independent measured datasets): {', '.join(high) or '—'}",
        f"  Medium confidence (1 independent measured dataset): {', '.join(medium) or '—'}",
    ]
    # Interpretive labels per category for the two halves of the claim.
    cat_label = {
        "receptors": "Neuromodulatory / synaptic receptors (→ neuromodulatory influences)",
        "excitability": "Intrinsic excitability genes (→ excitability profile)",
    }
    if by_cat:
        for cat in sorted(by_cat):
            lines.append(f"  {cat_label.get(cat, cat)}:")
            for fam in sorted(by_cat[cat]):
                lines.append(f"    - {fam}: {', '.join(by_cat[cat][fam])}")
    else:
        lines.append("  Expressed gene families:")
        for fam in sorted(summary["families_expressed"]):
            lines.append(f"    - {fam}: {', '.join(summary['families_expressed'][fam])}")
    lines.append(
        "  Note: tiers are detection-based; imputed Allen MERFISH genes are "
        "supporting evidence only (not independent of Allen scRNA)."
    )
    return "\n".join(lines)


def plot_evidence_dotplot(
    evidence: pd.DataFrame,
    config: dict[str, Any],
    *,
    cell_type: str,
    region: str | None = None,
    save_path: Path | str | None = None,
    output_dir: Path | str | None = None,
    max_genes: int = 60,
) -> Path | None:
    """Dot plot for one (cell_type, region): genes × datasets, size=frac, color=mean.

    Imputed Allen MERFISH cells are drawn with a hollow ring.
    """
    import matplotlib.pyplot as plt

    sub = target_evidence(evidence, cell_type=cell_type, region=region)
    sub = sub[sub["confidence_tier"] != "none"]
    if sub.empty:
        warnings.warn(
            f"No expressed genes for {cell_type!r} / {region!r}; skipping dot plot.",
            UserWarning,
            stacklevel=2,
        )
        return None
    if len(sub) > max_genes:
        sub = sub.head(max_genes)

    datasets = [k for k in DATASET_SPECS if f"{k}_mean" in sub.columns]
    genes = sub["gene"].tolist()
    gene_idx = {g: i for i, g in enumerate(genes)}

    cmap = config.get("output", {}).get("heatmap_cmap", "viridis")
    dpi = config.get("output", {}).get("dpi", 150)
    fig, ax = plt.subplots(figsize=(1.6 + 1.5 * len(datasets), 1.0 + 0.28 * len(genes)))

    all_means = sub[[f"{k}_mean" for k in datasets]].to_numpy(dtype=float)
    finite = all_means[np.isfinite(all_means)]
    vmax = float(np.nanpercentile(finite, 99)) if finite.size else 1.0
    vmax = vmax or 1.0

    for j, key in enumerate(datasets):
        means = sub[f"{key}_mean"].to_numpy(dtype=float)
        fracs = sub[f"{key}_frac"].to_numpy(dtype=float)
        sources = sub[f"{key}_source"].to_numpy()
        for g, m, fr, src in zip(genes, means, fracs, sources):
            if not np.isfinite(m):
                continue
            size = 30.0 + 320.0 * (0.0 if not np.isfinite(fr) else fr)
            ax.scatter(
                j, gene_idx[g], s=size, c=[m], cmap=cmap, vmin=0, vmax=vmax,
                edgecolors="red" if src == "imputed" else "black",
                linewidths=1.1 if src == "imputed" else 0.4,
            )

    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels([DATASET_SPECS[k]["label"] for k in datasets], rotation=30, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_ylim(-1, len(genes))
    ax.invert_yaxis()
    region_lbl = region or "(pooled)"
    ax.set_title(
        f"{cell_type} × {region_lbl}\nsize = fraction expressing, colour = mean log2(CPM+1); "
        "red ring = Allen imputed",
        fontsize=9,
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=vmax))
    fig.colorbar(sm, ax=ax, label="mean log2(CPM+1)", fraction=0.04, pad=0.02)
    fig.tight_layout()

    if save_path is None:
        out_root = resolve_output_dir(output_dir=output_dir, cfg=config if output_dir is None else None)
        safe_ct = "".join(c if c.isalnum() else "-" for c in cell_type).strip("-")
        safe_rg = (region or "pooled").replace("/", "-")
        save_path = Path(out_root) / f"evidence_dotplot_{safe_ct}_{safe_rg}.png"
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return save_path
