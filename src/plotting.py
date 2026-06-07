"""Heatmap and spatial plotting for receptor expression."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from scipy import stats

from src.config import get_figures_dir

IMPUTED_GENE_MARKER = "*"

_PROJECTION_AXES = {
    "coronal": ("x_ccf", "y_ccf"),
    "sagittal": ("z_ccf", "y_ccf"),
    "axial": ("x_ccf", "z_ccf"),
}

_SECTION_PROJECTION_AXES = {
    "coronal": ("x_section", "y_section"),
    "sagittal": ("z_section", "y_section"),
    "axial": ("x_section", "z_section"),
}


def _figsize(config: dict[str, Any], key: str, default: tuple[float, float]) -> tuple[float, float]:
    val = config["output"].get(key, list(default))
    return (float(val[0]), float(val[1]))


def _receptor_family_spans(
    genes: list[str],
    genes_flat: dict[str, str],
) -> list[tuple[str, int, int]]:
    """Return (family, start_col, end_col) inclusive for genes in display order."""
    spans: list[tuple[str, int, int]] = []
    idx = 0
    while idx < len(genes):
        family = genes_flat.get(genes[idx], "unknown")
        start = idx
        while idx < len(genes) and genes_flat.get(genes[idx]) == family:
            idx += 1
        spans.append((family, start, idx - 1))
    return spans


def _combined_heatmap_gene_labels(
    genes: list[str],
    config: dict[str, Any],
) -> list[str]:
    """Column labels; imputed Allen genes get a trailing asterisk when configured."""
    sources: dict[str, str] = config.get("_allen_gene_sources") or {}
    if not sources:
        return genes
    return [format_gene_label(g, sources.get(g, "measured")) for g in genes]


def _plot_combined_receptor_heatmap(
    mat: pd.DataFrame,
    config: dict[str, Any],
    title: str,
    figsize: tuple[float, float],
    save_path: Path | str,
) -> None:
    """Heatmap with cell types as rows, genes as columns, and ligand-family group labels."""
    cmap = config["output"].get("heatmap_cmap", "viridis")
    dpi = config["output"].get("dpi", 150)
    genes = list(mat.columns)
    display_genes = _combined_heatmap_gene_labels(genes, config)
    spans = _receptor_family_spans(genes, config["_genes_flat"])

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 14], hspace=0.04)
    ax_groups = fig.add_subplot(gs[0])
    ax_hm = fig.add_subplot(gs[1])

    data = mat.astype(float)
    values = data.to_numpy()
    mask = ~np.isfinite(values) if not np.isfinite(values).all() else None
    sns.heatmap(
        data,
        ax=ax_hm,
        cmap=cmap,
        mask=mask,
        cbar_kws={"label": "log2(CPM+1)"},
        xticklabels=display_genes,
        yticklabels=True,
    )
    ax_hm.set_xlabel("Receptor")
    ax_hm.set_ylabel(config["cell_type_level"])
    ax_hm.tick_params(axis="x", labelrotation=90, labelsize=7)
    ax_hm.tick_params(axis="y", labelsize=7)

    ax_groups.set_xlim(0, len(genes))
    ax_groups.set_ylim(0, 1)
    ax_groups.axis("off")
    for family, start, end in spans:
        x_left, x_right = start, end + 1
        cx = (x_left + x_right) / 2
        ax_groups.text(
            cx,
            0.55,
            family,
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
        )
        ax_groups.plot(
            [x_left, x_left, x_right, x_right],
            [0.05, 0.25, 0.25, 0.05],
            color="0.25",
            lw=1.0,
            clip_on=False,
        )

    fig.suptitle(title, y=0.99)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    agg_df: pd.DataFrame,
    title: str,
    config: dict[str, Any],
    save_path: Path | str | None = None,
    base_dir: Path | None = None,
    *,
    row_cluster: bool = True,
    col_cluster: bool = True,
    figsize: tuple[float, float] | None = None,
) -> None:
    """Seaborn clustermap wrapper; falls back to plain heatmap if too small."""
    if agg_df.empty:
        warnings.warn(f"Skipping empty heatmap: {title}", UserWarning, stacklevel=2)
        return

    cmap = config["output"].get("heatmap_cmap", "viridis")
    dpi = config["output"].get("dpi", 150)
    if figsize is None:
        figsize = _figsize(config, "figsize_heatmap", (14, 8))

    data = agg_df.astype(float)
    n_rows, n_cols = data.shape
    values = data.to_numpy()
    has_non_finite = not np.isfinite(values).all()

    can_cluster = n_rows >= 2 and n_cols >= 2 and not has_non_finite
    if has_non_finite and n_rows >= 2 and n_cols >= 2:
        warnings.warn(
            "Heatmap contains missing values (cell types absent in some brain areas); "
            "using plain heatmap without clustering.",
            UserWarning,
            stacklevel=2,
        )

    use_clustermap = can_cluster and (row_cluster or col_cluster)
    if use_clustermap:
        g = sns.clustermap(
            data,
            cmap=cmap,
            figsize=figsize,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            dendrogram_ratio=0.12,
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
        )
        g.fig.suptitle(title, y=1.02)
        fig = g.fig
    else:
        fig, ax = plt.subplots(figsize=figsize)
        mask = ~np.isfinite(values) if has_non_finite else None
        sns.heatmap(data, cmap=cmap, ax=ax, mask=mask)
        ax.set_title(title)

    if save_path:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def _scrna_heatmap_subtitle(config: dict[str, Any]) -> str:
    pools: dict[str, list[str]] = config.get("_scrna_pools") or {}
    if not pools:
        return ""
    parts = [
        f"{dissection} pools {', '.join(areas)}"
        for dissection, areas in pools.items()
    ]
    return f" [{'; '.join(parts)} — not distinguishable in scRNA]"


def format_gene_label(gene: str, source: str = "measured") -> str:
    """Append asterisk for imputed (non-measured) Allen MERFISH genes."""
    if source == "imputed":
        return f"{gene}{IMPUTED_GENE_MARKER}"
    return gene


def _modality_heatmap_subtitle(config: dict[str, Any]) -> str:
    mod = config.get("_dataset_modality")
    if mod == "zhuang":
        n = len(config.get("_zhuang_replicates_used") or [])
        return f" [Zhuang MERFISH, {n} replicates (mean), CCF parcellation]"
    if mod == "merfish":
        return " [Allen MERFISH, CCF parcellation]"
    return ""


def _merfish_heatmap_subtitle(config: dict[str, Any]) -> str:
    return _modality_heatmap_subtitle(config)


def plot_family_heatmap(
    family: str,
    gene_matrix: pd.DataFrame,
    config: dict[str, Any],
    base_dir: Path | None = None,
    output_dir: Path | str | None = None,
) -> Path:
    """Save heatmap_{family}.png under the configured figures directory."""
    figures_dir = get_figures_dir(config, base_dir, output_dir)
    out = figures_dir / f"heatmap_{family}.png"
    plot_heatmap(
        gene_matrix,
        title=(
            f"{family} receptors — mean log2(CPM+1)"
            f"{_merfish_heatmap_subtitle(config)}"
            f"{_scrna_heatmap_subtitle(config)}"
        ),
        config=config,
        save_path=out,
        base_dir=base_dir,
    )
    return out


def _cell_type_filter_subtitle(config: dict[str, Any]) -> str:
    patterns: list[str] = config.get("cell_type_name_filter") or []
    if not patterns:
        return ""
    return f" [name contains: {', '.join(patterns)}]"


def plot_combined_heatmap(
    all_genes_matrix: pd.DataFrame,
    config: dict[str, Any],
    base_dir: Path | None = None,
    output_dir: Path | str | None = None,
) -> Path:
    """Save heatmap_combined.png (cell types × receptor genes)."""
    figures_dir = get_figures_dir(config, base_dir, output_dir)
    out = figures_dir / "heatmap_combined.png"
    mat = all_genes_matrix.sort_index()
    n_cell_types, n_genes = mat.shape
    base = _figsize(config, "figsize_heatmap", (14, 8))
    figsize = (max(base[0], n_genes * 0.22), max(base[1], n_cell_types * 0.12))
    level = config["cell_type_level"]
    _plot_combined_receptor_heatmap(
        mat,
        config,
        title=(
            f"All {level}s × receptors — mean log2(CPM+1)"
            f"{_merfish_heatmap_subtitle(config)}"
            f"{_cell_type_filter_subtitle(config)}"
            f"{_scrna_heatmap_subtitle(config)}"
        ),
        figsize=figsize,
        save_path=out,
    )
    return out


def plot_spatial(
    coords_df: pd.DataFrame,
    expression: pd.Series,
    gene: str,
    projection: str,
    config: dict[str, Any],
    source_label: str = "",
    coord_prefix: str = "ccf",
    save_path: Path | str | None = None,
    base_dir: Path | None = None,
) -> None:
    """Scatter plot coloured by expression (CCF or section coordinates)."""
    axes_map = _PROJECTION_AXES if coord_prefix == "ccf" else _SECTION_PROJECTION_AXES
    if projection not in axes_map:
        raise ValueError(f"projection must be one of {list(axes_map)}")

    ax1, ax2 = axes_map[projection]
    common = coords_df.index.intersection(expression.index)
    coords = coords_df.loc[common]
    expr = expression.loc[common].astype(float)

    cmap = config["output"].get("spatial_cmap", "magma")
    figsize = _figsize(config, "figsize_spatial", (10, 10))
    dpi = config["output"].get("dpi", 150)

    fig, ax = plt.subplots(figsize=figsize)
    positive = expr[expr > 0]
    vmax = float(np.percentile(positive, 99)) if len(positive) > 0 else float(expr.max() or 1)

    sc = ax.scatter(
        coords[ax1],
        coords[ax2],
        c=expr,
        cmap=cmap,
        s=0.3,
        alpha=0.6,
        vmin=0,
        vmax=vmax,
        rasterized=True,
    )
    plt.colorbar(sc, ax=ax, label="log2(CPM+1)")

    title = f"{gene} — MERFISH {projection}"
    if source_label == "imputed":
        title += " [imputed]"
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.invert_yaxis()

    if save_path:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def plot_family_spatial_panel(
    family_results: dict[str, dict[str, Any]],
    family: str,
    coords_df: pd.DataFrame,
    config: dict[str, Any],
    base_dir: Path | None = None,
    output_dir: Path | str | None = None,
) -> Path | None:
    """
    Multi-panel grid: one row per gene, one column per projection.

    family_results: {gene: {'expression': Series, 'source': str}}
    """
    genes = list(family_results.keys())
    if not genes:
        return None

    projections = ["coronal", "sagittal", "axial"]
    nrows, ncols = len(genes), len(projections)
    figsize = _figsize(config, "figsize_spatial", (10, 10))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize[0] * ncols * 0.35, figsize[1] * nrows * 0.35),
        squeeze=False,
    )
    cmap = config["output"].get("spatial_cmap", "magma")
    dpi = config["output"].get("dpi", 150)

    for i, gene in enumerate(genes):
        info = family_results[gene]
        expr = info["expression"]
        source = info.get("source", "")
        for j, proj in enumerate(projections):
            ax = axes[i, j]
            ax1, ax2 = _PROJECTION_AXES[proj]
            common = coords_df.index.intersection(expr.index)
            coords = coords_df.loc[common]
            e = expr.loc[common].astype(float)
            positive = e[e > 0]
            vmax = float(np.percentile(positive, 99)) if len(positive) > 0 else float(e.max() or 1)
            ax.scatter(
                coords[ax1],
                coords[ax2],
                c=e,
                cmap=cmap,
                s=0.2,
                alpha=0.5,
                vmin=0,
                vmax=vmax,
                rasterized=True,
            )
            subtitle = f"{gene} {proj}"
            if source == "imputed":
                subtitle += " [imp]"
            ax.set_title(subtitle, fontsize=8)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle(f"{family} — spatial expression (CCF)", fontsize=12)
    fig.tight_layout()

    figures_dir = get_figures_dir(config, base_dir, output_dir)
    out = figures_dir / f"spatial_panel_{family}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_crossref_scatter(
    merged: pd.DataFrame,
    config: dict[str, Any],
    *,
    title: str,
    save_path: Path | str,
    color_by: str = "brain_area",
) -> None:
    """
    Scatter Allen MERFISH vs Zhuang mean expression.

    Imputed Allen genes are drawn with hollow markers; measured use filled circles.
    """
    if merged.empty:
        warnings.warn(f"Skipping empty cross-reference plot: {title}", UserWarning, stacklevel=2)
        return

    dpi = config["output"].get("dpi", 150)
    figsize = _figsize(config, "figsize_heatmap", (14, 8))

    fig, ax = plt.subplots(figsize=figsize)
    areas = sorted(merged[color_by].dropna().unique())
    palette = sns.color_palette("tab10", n_colors=max(len(areas), 1))
    area_colors = dict(zip(areas, palette))

    for area in areas:
        sub = merged[merged[color_by] == area]
        for source, marker, facecolors in (
            ("measured", "o", "auto"),
            ("imputed", "o", "none"),
        ):
            pts = sub[sub["allen_source"] == source]
            if pts.empty:
                continue
            ax.scatter(
                pts["allen_expression"],
                pts["zhuang_expression"],
                c=[area_colors[area]],
                label=f"{area}" if source == "measured" else None,
                s=28 if source == "measured" else 40,
                alpha=0.75,
                marker=marker,
                edgecolors=[area_colors[area]] if source == "imputed" else "none",
                linewidths=0.8 if source == "imputed" else 0,
            )

    pearson_r, pearson_p = stats.pearsonr(
        merged["allen_expression"],
        merged["zhuang_expression"],
    )
    spearman_r, _ = stats.spearmanr(
        merged["allen_expression"],
        merged["zhuang_expression"],
    )

    lim_lo = float(
        min(merged["allen_expression"].min(), merged["zhuang_expression"].min()),
    )
    lim_hi = float(
        max(merged["allen_expression"].max(), merged["zhuang_expression"].max()),
    )
    pad = (lim_hi - lim_lo) * 0.05 or 0.1
    lims = (lim_lo - pad, lim_hi + pad)
    ax.plot(lims, lims, "k--", alpha=0.35, linewidth=1, zorder=0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    n_imputed = int((merged["allen_source"] == "imputed").sum())
    imputed_note = (
        f"; hollow = imputed Allen ({n_imputed} points)"
        if n_imputed
        else ""
    )
    ax.set_xlabel("Allen MERFISH — mean log2(CPM+1)")
    ax.set_ylabel("Zhuang MERFISH — mean log2(CPM+1)")
    ax.set_title(
        f"{title}\n"
        f"n={len(merged)} cell_type×region×gene; "
        f"r={pearson_r:.3f} (p={pearson_p:.2g}), "
        f"ρ={spearman_r:.3f}{imputed_note}",
        fontsize=11,
    )
    ax.set_aspect("equal", adjustable="box")
    if areas:
        ax.legend(title=color_by, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)

    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_crossref_family_scatters(
    merged: pd.DataFrame,
    config: dict[str, Any],
    output_dir: Path | str | None = None,
) -> list[Path]:
    """Save overall and per-family Allen vs Zhuang correlation scatter plots."""
    figures_dir = get_figures_dir(config, output_dir=output_dir)
    saved: list[Path] = []

    overall = figures_dir / "crossref_allen_zhuang.png"
    plot_crossref_scatter(
        merged,
        config,
        title="Allen MERFISH vs Zhuang MERFISH",
        save_path=overall,
    )
    saved.append(overall)

    for family in config["_families"]:
        sub = merged[merged["family"] == family]
        if sub.empty:
            warnings.warn(
                f"No cross-reference data for family {family!r}; skipping scatter.",
                UserWarning,
                stacklevel=2,
            )
            continue
        out = figures_dir / f"crossref_allen_zhuang_{family}.png"
        plot_crossref_scatter(
            sub,
            config,
            title=f"Allen vs Zhuang — {family}",
            save_path=out,
        )
        saved.append(out)

    return saved
