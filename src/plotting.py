"""Heatmap and spatial plotting for receptor expression."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from src.config import get_figures_dir

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


def plot_heatmap(
    agg_df: pd.DataFrame,
    title: str,
    config: dict[str, Any],
    save_path: Path | str | None = None,
    base_dir: Path | None = None,
) -> None:
    """Seaborn clustermap wrapper; falls back to plain heatmap if too small."""
    if agg_df.empty:
        warnings.warn(f"Skipping empty heatmap: {title}", UserWarning, stacklevel=2)
        return

    cmap = config["output"].get("heatmap_cmap", "viridis")
    dpi = config["output"].get("dpi", 150)
    figsize = _figsize(config, "figsize_heatmap", (14, 8))

    data = agg_df.astype(float)
    n_rows, n_cols = data.shape

    if n_rows >= 2 and n_cols >= 2:
        g = sns.clustermap(
            data,
            cmap=cmap,
            figsize=figsize,
            row_cluster=True,
            col_cluster=True,
            dendrogram_ratio=0.12,
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
        )
        g.fig.suptitle(title, y=1.02)
        fig = g.fig
    else:
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(data, cmap=cmap, ax=ax)
        ax.set_title(title)

    if save_path:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def plot_family_heatmap(
    family: str,
    gene_matrix: pd.DataFrame,
    config: dict[str, Any],
    base_dir: Path | None = None,
) -> Path:
    """Save figures/heatmap_{family}.png."""
    figures_dir = get_figures_dir(config, base_dir)
    out = figures_dir / f"heatmap_{family}.png"
    plot_heatmap(
        gene_matrix,
        title=f"{family} receptors — mean log2(CPM+1)",
        config=config,
        save_path=out,
        base_dir=base_dir,
    )
    return out


def plot_combined_heatmap(
    all_genes_matrix: pd.DataFrame,
    config: dict[str, Any],
    base_dir: Path | None = None,
) -> Path:
    """Save figures/heatmap_combined.png (genes x top cell types)."""
    figures_dir = get_figures_dir(config, base_dir)
    out = figures_dir / "heatmap_combined.png"
    # Rows = cell types, cols = genes — transpose for clustermap readability
    mat = all_genes_matrix
    plot_heatmap(
        mat,
        title="All receptors — top variable cell types",
        config=config,
        save_path=out,
        base_dir=base_dir,
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

    figures_dir = get_figures_dir(config, base_dir)
    out = figures_dir / f"spatial_panel_{family}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out
