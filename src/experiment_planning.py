"""NE/ACh experiment planning from synthesis evidence tables.

Loads prior synthesis outputs, ranks pharmacology targets, predicts combined
transmitter effects at low vs high concentration, writes figures and a text plan.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

from src.functional_analysis import infer_coarse_type, pick_expression_value
from src.pharmacology import (
    TIER_ORDER,
    get_pharmacology,
    coupling_sign,
    transmitter_weight,
)

EXPRESSED_TIERS = frozenset({"high", "medium"})
FAMILY_LABEL = {"noradrenaline": "NE", "acetylcholine": "ACh"}
TRANSMITTER_BY_FAMILY = {"noradrenaline": "NE", "acetylcholine": "ACh"}
COUPLING_COLORS = {"Gi": "#4C72B0", "Gq": "#C44E52", "Gs": "#55A868", "ionotropic": "#8172B2"}
COUPLING_LEGEND = (
    "Bar colours = G-protein / channel coupling: "
    "Gi (blue), Gq (red), Gs (green), ionotropic (purple)"
)
AXIS_LABELS = ["excitability", "rin", "ih", "mcurrent", "adaptation"]


@dataclass
class PlanningSpec:
    """Scope and presentation for an experiment-planning report."""

    brain_areas: tuple[str, ...] = ("VISp", "V2M")
    cell_type_patterns: tuple[str, ...] | None = None
    coarse_types: tuple[str, ...] | None = None
    per_supertype: bool = False
    family_order: tuple[str, ...] = ("acetylcholine", "noradrenaline")
    compare_groups: tuple[str, ...] | None = None
    report_title: str = "NE & ACh RECEPTOR EXPRESSION — EXPERIMENT PLAN"
    plan_filename: str = "EXPERIMENT_PLAN_NE_ACh.txt"
    region_label: str = "Regions: VISp (V1) vs V2M (VISpm + VISam + RSPagl unweighted rollup)"
    include_region_delta: bool = True
    include_msc_handover: bool = False
    msc_project_title: str = ""
    extra_caveats: tuple[str, ...] = field(default_factory=tuple)


DEFAULT_VISp_V2M_SPEC = PlanningSpec(
    family_order=("noradrenaline", "acetylcholine"),
    compare_groups=None,
    include_region_delta=True,
)

V2M_L5_IT_ET_SPEC = PlanningSpec(
    brain_areas=("V2M",),
    cell_type_patterns=("L5 IT CTX", "L5 ET CTX"),
    coarse_types=("L5 IT", "L5 ET"),
    per_supertype=True,
    family_order=("acetylcholine", "noradrenaline"),
    compare_groups=("L5 IT", "L5 ET"),
    report_title="V2M L5 IT & L5 ET — ACh / NE EXPERIMENT PLAN (MSc project)",
    plan_filename="EXPERIMENT_PLAN_V2M_L5_ACh_NE.txt",
    region_label="Region: V2M only (VISpm + VISam + RSPagl unweighted rollup)",
    include_region_delta=False,
    include_msc_handover=True,
    msc_project_title="V2M L5 pyramidal neuromodulation — ACh (primary) and NE (secondary)",
    extra_caveats=(
        "V2M is a synthetic rollup; confirm recorded cells lie in VISpm, VISam, or RSPagl.",
        "L5 IT vs L5 ET are Allen supertypes — validate by morphology/projection during recording.",
        "Thalamic (LP/ATN) afferents to V2M L5 ET tuft are nicotinic-facilitated (see receptor_excitability.md §3).",
    ),
)


def find_latest_synthesis_run(
    exploration_root: Path | str,
    cell_type_level: str,
) -> Path:
    """Most recent ``{timestamp}_{level}_synthesis`` folder with an evidence table."""
    root = Path(exploration_root).expanduser()
    suffix = f"_{cell_type_level}_synthesis"
    candidates = [
        d for d in root.iterdir()
        if d.is_dir() and d.name.endswith(suffix)
        and ((d / "evidence_table.parquet").is_file() or (d / "evidence_table.csv").is_file())
    ]
    if not candidates:
        raise FileNotFoundError(
            f"No synthesis run found under {root} for level {cell_type_level!r}. "
            "Run notebook 05_synthesis first."
        )
    return sorted(candidates, key=lambda p: p.name, reverse=True)[0]


def load_evidence_table(run_dir: Path | str) -> pd.DataFrame:
    run_dir = Path(run_dir)
    pq = run_dir / "evidence_table.parquet"
    csv = run_dir / "evidence_table.csv"
    if pq.is_file():
        return pd.read_parquet(pq)
    if csv.is_file():
        return pd.read_csv(csv, low_memory=False)
    raise FileNotFoundError(f"No evidence table in {run_dir}")


def filter_planning_evidence(
    evidence: pd.DataFrame,
    config: dict[str, Any],
    *,
    spec: PlanningSpec | None = None,
    families: tuple[str, ...] | None = None,
    brain_areas: tuple[str, ...] | None = None,
    min_tier: str = "medium",
) -> pd.DataFrame:
    """High/medium NE or ACh rows for selected regions and pyramidal cell types."""
    spec = spec or DEFAULT_VISp_V2M_SPEC
    families = families or spec.family_order
    brain_areas = brain_areas or spec.brain_areas
    allowed_tiers = {"high", "medium"} if min_tier == "medium" else {"high"}
    sub = evidence[
        evidence["family"].isin(families)
        & evidence["brain_area"].isin(brain_areas)
        & evidence["confidence_tier"].isin(allowed_tiers)
    ].copy()

    patterns = spec.cell_type_patterns
    if patterns:
        mask = pd.Series(False, index=sub.index)
        for token in patterns:
            mask |= sub["cell_type"].str.contains(token, na=False, regex=False)
        sub = sub[mask]
    else:
        name_filter = config.get("cell_type_name_filter") or []
        if name_filter:
            mask = pd.Series(False, index=sub.index)
            for token in name_filter:
                mask |= sub["cell_type"].str.contains(token, na=False, regex=False)
            sub = sub[mask]

    sub["coarse_type"] = sub["cell_type"].map(lambda s: infer_coarse_type(str(s)))
    sub = sub[sub["coarse_type"].notna()].copy()
    if spec.coarse_types:
        sub = sub[sub["coarse_type"].isin(spec.coarse_types)].copy()

    expr_vals: list[float] = []
    expr_ds: list[str | None] = []
    for _, row in sub.iterrows():
        val, ds, _src = pick_expression_value(row, config)
        expr_vals.append(val)
        expr_ds.append(ds)
    sub["expression"] = expr_vals
    sub["expression_dataset"] = expr_ds
    return sub


def aggregate_expression_summary(
    sub: pd.DataFrame,
    *,
    per_supertype: bool = False,
) -> pd.DataFrame:
    """Aggregate to coarse_type × region × gene, or keep per supertype."""
    if per_supertype:
        return _aggregate_supertype_summary(sub)
    return aggregate_coarse_summary(sub)


def _aggregate_supertype_summary(sub: pd.DataFrame) -> pd.DataFrame:
    """One row per Allen supertype × brain_area × gene."""
    rows: list[dict[str, Any]] = []
    group_cols = ["cell_type", "coarse_type", "brain_area", "gene"]
    for key, chunk in sub.groupby(group_cols, observed=True):
        ct, coarse, area, gene = key
        chunk = chunk.assign(_tier_ord=chunk["confidence_tier"].map(TIER_ORDER))
        chunk = chunk.sort_values(["_tier_ord", "expression"], ascending=[True, False])
        best = chunk.iloc[0]
        rows.append({
            "cell_type": ct,
            "coarse_type": coarse,
            "brain_area": area,
            "gene": gene,
            "family": best["family"],
            "confidence_tier": best["confidence_tier"],
            "expression": float(chunk["expression"].max()),
            "expression_dataset": best["expression_dataset"],
            "top_cell_type": ct,
            "n_supertypes": 1,
        })
    return pd.DataFrame(rows)


def aggregate_coarse_summary(sub: pd.DataFrame) -> pd.DataFrame:
    """One row per coarse_type × brain_area × gene (best tier, max expression)."""
    rows: list[dict[str, Any]] = []
    for (coarse, area, gene), chunk in sub.groupby(
        ["coarse_type", "brain_area", "gene"], observed=True,
    ):
        chunk = chunk.assign(_tier_ord=chunk["confidence_tier"].map(TIER_ORDER))
        chunk = chunk.sort_values(["_tier_ord", "expression"], ascending=[True, False])
        best = chunk.iloc[0]
        top_ct = chunk.loc[chunk["expression"].idxmax(), "cell_type"]
        rows.append({
            "coarse_type": coarse,
            "brain_area": area,
            "gene": gene,
            "family": best["family"],
            "confidence_tier": best["confidence_tier"],
            "expression": float(chunk["expression"].max()),
            "expression_dataset": chunk.loc[chunk["expression"].idxmax(), "expression_dataset"],
            "top_cell_type": top_ct,
            "n_supertypes": int(chunk["cell_type"].nunique()),
        })
    return pd.DataFrame(rows)


def score_pharmacology_priority(row: pd.Series) -> float:
    tier_bonus = {"high": 3.0, "medium": 1.5}.get(row["confidence_tier"], 0.0)
    expr = float(row["expression"]) if pd.notna(row["expression"]) else 0.0
    meta = get_pharmacology(str(row["gene"]))
    coupling_bonus = {"Gq": 0.5, "Gs": 0.4, "Gi": 0.3, "ionotropic": 0.35}.get(
        meta.coupling if meta else "", 0.0,
    )
    return tier_bonus + 0.25 * expr + coupling_bonus


def build_receptor_order_table(summary: pd.DataFrame) -> pd.DataFrame:
    """Ranked pharmacology candidates with compound suggestions."""
    rows: list[dict[str, Any]] = []
    for _, row in summary.iterrows():
        meta = get_pharmacology(str(row["gene"]))
        if meta is None:
            continue
        rows.append({
            "cell_type": row.get("cell_type", row.get("top_cell_type")),
            "coarse_type": row["coarse_type"],
            "brain_area": row["brain_area"],
            "gene": row["gene"],
            "receptor": meta.receptor,
            "coupling": meta.coupling,
            "confidence_tier": row["confidence_tier"],
            "expression": row["expression"],
            "top_cell_type": row["top_cell_type"],
            "priority_score": score_pharmacology_priority(row),
            "intrinsic_effect": meta.intrinsic_effect,
            "rin": meta.rin_effect,
            "ih": meta.ih_effect,
            "mcurrent": meta.mcurrent_effect,
            "firing": meta.firing_effect,
            "agonists": "; ".join(meta.agonists),
            "antagonists": "; ".join(meta.antagonists),
            "notes": meta.notes,
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    sort_cols = ["coarse_type", "brain_area", "priority_score"]
    if "cell_type" in out.columns and out["cell_type"].nunique() > out["coarse_type"].nunique():
        sort_cols = ["coarse_type", "cell_type", "priority_score"]
    return out.sort_values(sort_cols, ascending=[True, True, False])


def compute_transmitter_scenario(
    summary: pd.DataFrame,
    *,
    coarse_type: str,
    brain_area: str,
    family: str,
    concentration: Literal["low", "high"],
) -> dict[str, Any]:
    """Weighted coupling activity for combined NE or ACh at low/high tone."""
    transmitter = TRANSMITTER_BY_FAMILY[family]
    sub = summary[
        (summary["coarse_type"] == coarse_type)
        & (summary["brain_area"] == brain_area)
        & (summary["family"] == family)
    ]
    coupling_totals: dict[str, float] = {c: 0.0 for c in COUPLING_COLORS}
    axis_totals = {a: 0.0 for a in AXIS_LABELS}
    active_genes: list[str] = []

    for _, row in sub.iterrows():
        gene = str(row["gene"])
        meta = get_pharmacology(gene)
        if meta is None:
            continue
        aff = transmitter_weight(gene, transmitter, concentration)
        if aff <= 0:
            continue
        expr_w = float(row["expression"]) if pd.notna(row["expression"]) else 1.0
        tier_w = 1.0 if row["confidence_tier"] == "high" else 0.65
        weight = aff * tier_w * (1.0 + 0.1 * expr_w)
        coupling_totals[meta.coupling] += weight
        signs = coupling_sign(meta.coupling)
        for axis in AXIS_LABELS:
            axis_totals[axis] += weight * signs[axis]
        active_genes.append(gene)

    total = sum(coupling_totals.values()) or 1.0
    coupling_frac = {k: v / total for k, v in coupling_totals.items()}
    return {
        "coarse_type": coarse_type,
        "brain_area": brain_area,
        "family": family,
        "transmitter": transmitter,
        "concentration": concentration,
        "coupling_totals": coupling_totals,
        "coupling_fraction": coupling_frac,
        "axis_scores": axis_totals,
        "active_genes": sorted(set(active_genes)),
    }


def build_scenario_table(
    summary: pd.DataFrame,
    coarse_types: list[str],
    brain_areas: list[str],
    families: tuple[str, ...] = ("noradrenaline", "acetylcholine"),
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for family in families:
        for coarse in coarse_types:
            for area in brain_areas:
                for conc in ("low", "high"):
                    scen = compute_transmitter_scenario(
                        summary, coarse_type=coarse, brain_area=area,
                        family=family, concentration=conc,
                    )
                    rec = {
                        "family": family,
                        "transmitter": scen["transmitter"],
                        "coarse_type": coarse,
                        "brain_area": area,
                        "concentration": conc,
                        "active_genes": ", ".join(scen["active_genes"]),
                        **{f"coupling_{k}": v for k, v in scen["coupling_fraction"].items()},
                        **{f"axis_{k}": v for k, v in scen["axis_scores"].items()},
                    }
                    rows.append(rec)
    return pd.DataFrame(rows)


def coarse_summary_for_scenarios(summary: pd.DataFrame) -> pd.DataFrame:
    """Collapse supertype rows to coarse_type × brain_area × gene (max expression)."""
    rows: list[dict[str, Any]] = []
    for (coarse, area, gene), chunk in summary.groupby(
        ["coarse_type", "brain_area", "gene"], observed=True,
    ):
        chunk = chunk.assign(_tier_ord=chunk["confidence_tier"].map(TIER_ORDER))
        chunk = chunk.sort_values(["_tier_ord", "expression"], ascending=[True, False])
        best = chunk.iloc[0]
        top = chunk.loc[chunk["expression"].idxmax()]
        rows.append({
            "coarse_type": coarse,
            "brain_area": area,
            "gene": gene,
            "family": best["family"],
            "confidence_tier": best["confidence_tier"],
            "expression": float(chunk["expression"].max()),
            "expression_dataset": top["expression_dataset"],
            "top_cell_type": top.get("cell_type", top.get("top_cell_type")),
            "n_supertypes": int(chunk["cell_type"].nunique()) if "cell_type" in chunk else 1,
        })
    return pd.DataFrame(rows)


def _compare_group_delta(
    summary: pd.DataFrame,
    gene: str,
    group_a: str,
    group_b: str,
    *,
    brain_area: str | None = None,
) -> float:
    """Mean expression in group_b − group_a (e.g. L5 ET − L5 IT)."""
    a = summary[(summary["gene"] == gene) & (summary["coarse_type"] == group_a)]
    b = summary[(summary["gene"] == gene) & (summary["coarse_type"] == group_b)]
    if brain_area:
        a = a[a["brain_area"] == brain_area]
        b = b[b["brain_area"] == brain_area]
    if a.empty or b.empty:
        return np.nan
    return float(b["expression"].mean() - a["expression"].mean())


def _region_delta(summary: pd.DataFrame, gene: str, coarse: str) -> float:
    visp = summary[
        (summary["gene"] == gene) & (summary["coarse_type"] == coarse)
        & (summary["brain_area"] == "VISp")
    ]["expression"]
    v2m = summary[
        (summary["gene"] == gene) & (summary["coarse_type"] == coarse)
        & (summary["brain_area"] == "V2M")
    ]["expression"]
    if visp.empty or v2m.empty:
        return np.nan
    return float(v2m.max() - visp.max())


def plot_expression_heatmap(
    summary: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
    spec: PlanningSpec | None = None,
) -> Path:
    """Expression heatmap: coarse types or individual supertypes × genes."""
    spec = spec or DEFAULT_VISp_V2M_SPEC
    sub = summary[summary["family"] == family].copy()
    if sub.empty:
        raise ValueError(f"No summary rows for family {family!r}")

    genes = sorted(sub["gene"].unique(), key=lambda g: (
        -sub.loc[sub["gene"] == g, "expression"].max(),
        g,
    ))

    if spec.per_supertype and "cell_type" in sub.columns:
        row_order = sorted(sub["cell_type"].unique())
        mat = pd.DataFrame(index=row_order, columns=genes, dtype=float)
        for _, row in sub.iterrows():
            mat.loc[row["cell_type"], row["gene"]] = row["expression"]
        row_labels = list(row_order)
        title_suffix = f"V2M L5 supertypes — {FAMILY_LABEL.get(family, family)}"
    else:
        coarse_order = list(spec.coarse_types or ["L2/3 IT", "L5 IT", "L5 ET"])
        areas = list(spec.brain_areas)
        index = pd.MultiIndex.from_product([coarse_order, areas], names=["coarse_type", "brain_area"])
        mat = pd.DataFrame(index=index, columns=genes, dtype=float)
        for _, row in sub.iterrows():
            key = (row["coarse_type"], row["brain_area"])
            if key not in mat.index:
                continue
            mat.loc[key, row["gene"]] = row["expression"]
        row_labels = [f"{c} · {a}" for c, a in mat.index]
        title_suffix = f"{FAMILY_LABEL.get(family, family)} — coarse types"

    fig_h = max(4.5, 0.32 * len(mat))
    fig_w = max(10, 0.45 * len(genes) + 4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = config["output"].get("heatmap_cmap", "viridis")
    sns.heatmap(
        mat.astype(float),
        ax=ax,
        cmap=cmap,
        cbar_kws={"label": "log2(CPM+1)"},
        linewidths=0.5,
        linecolor="white",
        yticklabels=row_labels,
        xticklabels=genes,
    )
    ax.set_title(title_suffix)
    ax.set_xlabel("Gene")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", labelsize=8)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(fig)
    return path


def _add_coupling_legend(ax: plt.Axes) -> None:
    handles = [
        mpatches.Patch(color=c, label=k) for k, c in COUPLING_COLORS.items()
    ]
    ax.legend(handles=handles, title="Coupling", loc="lower right", fontsize=8)


def plot_priority_bars(
    order_table: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
    top_n: int = 12,
    title_suffix: str = "",
) -> Path:
    sub = order_table[
        order_table["gene"].map(
            lambda g: (get_pharmacology(g) or None) is not None
            and get_pharmacology(g).family == family  # type: ignore[union-attr]
        )
    ]
    if sub.empty:
        raise ValueError(f"No order rows for {family}")

    plot_df = sub.nlargest(top_n, "priority_score").copy()
    if "cell_type" in plot_df.columns and plot_df["cell_type"].nunique() > 1:
        plot_df["label"] = plot_df.apply(
            lambda r: f"{r['gene']} ({r['cell_type']})", axis=1,
        )
    else:
        plot_df["label"] = plot_df.apply(
            lambda r: f"{r['gene']} ({r['coarse_type']}, {r['brain_area']})", axis=1,
        )

    fig, ax = plt.subplots(figsize=(10, max(5, 0.35 * len(plot_df))))
    colors = plot_df["coupling"].map(COUPLING_COLORS).fillna("#999999")
    ax.barh(plot_df["label"], plot_df["priority_score"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("Priority score (tier + expression + coupling)")
    fam_lbl = FAMILY_LABEL.get(family, family)
    ax.set_title(f"{fam_lbl} — pharmacology ordering{title_suffix}")
    _add_coupling_legend(ax)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(fig)
    return path


def plot_it_et_delta(
    summary: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
    spec: PlanningSpec,
) -> Path:
    """L5 ET − L5 IT mean expression at V2M (per gene)."""
    if not spec.compare_groups or len(spec.compare_groups) != 2:
        raise ValueError("spec.compare_groups must be (group_a, group_b)")
    group_a, group_b = spec.compare_groups
    sub = summary[summary["family"] == family]
    genes = sorted(sub["gene"].unique())
    area = spec.brain_areas[0] if len(spec.brain_areas) == 1 else None
    rows: list[dict[str, Any]] = []
    for gene in genes:
        d = _compare_group_delta(summary, gene, group_a, group_b, brain_area=area)
        if np.isfinite(d):
            rows.append({"gene": gene, "delta_et_minus_it": d})
    plot_df = pd.DataFrame(rows)
    if plot_df.empty:
        raise ValueError("No IT vs ET deltas to plot")

    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(genes)), 4.5))
    colors = plot_df["gene"].map(
        lambda g: COUPLING_COLORS.get(get_pharmacology(g).coupling if get_pharmacology(g) else "", "#999")
    )
    ax.bar(plot_df["gene"], plot_df["delta_et_minus_it"], color=colors)
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xlabel("Gene")
    ax.set_ylabel("L5 ET − L5 IT expression (log2 CPM+1)")
    fam_lbl = FAMILY_LABEL.get(family, family)
    ax.set_title(f"{fam_lbl} — V2M subtype difference (positive = higher in L5 ET)")
    ax.tick_params(axis="x", rotation=45)
    _add_coupling_legend(ax)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(fig)
    return path


def plot_subtype_dot_comparison(
    summary: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
) -> Path:
    """Dot plot: each supertype × gene expression, coloured by coarse IT/ET."""
    sub = summary[summary["family"] == family].copy()
    if sub.empty or "cell_type" not in sub.columns:
        raise ValueError("Need per-supertype summary for dot comparison")

    palette = {"L5 IT": "#4C72B0", "L5 ET": "#C44E52"}
    sub = sub.sort_values(["gene", "coarse_type", "cell_type"])
    fig_h = max(5, 0.25 * sub["cell_type"].nunique())
    fig, ax = plt.subplots(figsize=(max(9, 0.4 * sub["gene"].nunique() + 3), fig_h))
    sns.scatterplot(
        data=sub,
        x="gene",
        y="cell_type",
        size="expression",
        hue="coarse_type",
        palette=palette,
        sizes=(40, 400),
        ax=ax,
        legend="brief",
    )
    fam_lbl = FAMILY_LABEL.get(family, family)
    ax.set_title(f"{fam_lbl} — V2M L5 supertype expression (dot size ∝ log2 CPM+1)")
    ax.tick_params(axis="x", rotation=45)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(fig)
    return path


def plot_region_delta(
    summary: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
) -> Path:
    sub = summary[summary["family"] == family]
    coarse_order = ["L2/3 IT", "L5 IT", "L5 ET"]
    genes = sorted(sub["gene"].unique())
    rows: list[dict[str, Any]] = []
    for coarse in coarse_order:
        for gene in genes:
            d = _region_delta(sub, gene, coarse)
            if np.isfinite(d):
                rows.append({"coarse_type": coarse, "gene": gene, "delta_v2m_visp": d})
    plot_df = pd.DataFrame(rows)
    if plot_df.empty:
        raise ValueError("No region deltas to plot")

    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(genes)), 4.5))
    sns.barplot(
        data=plot_df,
        x="gene",
        y="delta_v2m_visp",
        hue="coarse_type",
        ax=ax,
        palette="Set2",
    )
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xlabel("Gene")
    ax.set_ylabel("V2M − VISp expression (log2 CPM+1)")
    fam_lbl = FAMILY_LABEL.get(family, family)
    ax.set_title(f"{fam_lbl} — regional difference (positive = higher in V2M rollup)")
    ax.tick_params(axis="x", rotation=45)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(fig)
    return path


def plot_coupling_scenarios(
    scenario_table: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
) -> Path:
    sub = scenario_table[scenario_table["family"] == family].copy()
    if sub.empty:
        raise ValueError(f"No scenarios for {family}")

    sub["panel"] = sub["coarse_type"] + " · " + sub["brain_area"]
    coupling_cols = [c for c in sub.columns if c.startswith("coupling_")]
    plot_rows: list[dict[str, Any]] = []
    for _, row in sub.iterrows():
        for col in coupling_cols:
            coupling = col.replace("coupling_", "")
            plot_rows.append({
                "panel": row["panel"],
                "concentration": row["concentration"],
                "coupling": coupling,
                "fraction": row[col],
            })
    plot_df = pd.DataFrame(plot_rows)
    g = sns.catplot(
        data=plot_df,
        kind="bar",
        x="panel",
        y="fraction",
        hue="coupling",
        col="concentration",
        palette=COUPLING_COLORS,
        height=4.5,
        aspect=1.4,
        legend=True,
    )
    fam_lbl = FAMILY_LABEL.get(family, family)
    g.fig.suptitle(
        f"{fam_lbl} combined receptor activation by coupling class (expression × affinity)",
        y=1.02,
    )
    g.set_xticklabels(rotation=45, ha="right")
    g.set_axis_labels("", "Weighted coupling fraction")
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    g.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(g.fig)
    return path


def plot_axis_scenarios(
    scenario_table: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
) -> Path:
    sub = scenario_table[scenario_table["family"] == family].copy()
    axis_cols = [c for c in sub.columns if c.startswith("axis_")]
    if sub.empty or not axis_cols:
        raise ValueError(f"No axis scenarios for {family}")

    sub["panel"] = sub["coarse_type"] + " · " + sub["brain_area"]
    mat = sub.set_index(["panel", "concentration"])[axis_cols]
    mat.columns = [c.replace("axis_", "") for c in mat.columns]
    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, 0.35 * len(mat) // 2)), sharey=True)
    for ax, conc in zip(axes, ["low", "high"]):
        chunk = mat.xs(conc, level="concentration")
        sns.heatmap(
            chunk.astype(float),
            ax=ax,
            cmap="RdBu_r",
            center=0,
            cbar=ax is axes[-1],
            yticklabels=True,
            xticklabels=True,
        )
        ax.set_title(f"{conc.upper()} {FAMILY_LABEL.get(family, family)} tone")
        ax.set_xlabel("Predicted intrinsic axis")
    fig.suptitle("Combined transmitter effect on intrinsic excitability axes (+ = increase / excitation)")
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(fig)
    return path


def _interpret_axis_scores(scores: dict[str, float]) -> list[str]:
    lines: list[str] = []
    if scores.get("excitability", 0) > 0.5:
        lines.append("Net shift toward higher intrinsic excitability / depolarisation.")
    elif scores.get("excitability", 0) < -0.3:
        lines.append("Net inhibitory/shunting influence on intrinsic excitability.")
    if scores.get("rin", 0) > 0.4:
        lines.append("Expect increased input resistance (α2-like or Gq ↓Kir2/M-current).")
    elif scores.get("rin", 0) < -0.3:
        lines.append("Expect decreased input resistance (Gs/β-HCN or ionotropic shunt).")
    if scores.get("ih", 0) > 0.4:
        lines.append("HCN / I_h axis biased toward increase (β1/Gs rightward shift).")
    elif scores.get("ih", 0) < -0.4:
        lines.append("HCN / I_h axis biased toward decrease (α2 PLC-PKC leftward shift).")
    if scores.get("mcurrent", 0) < -0.4:
        lines.append("M-current suppression expected → loss of spike-frequency adaptation.")
    if scores.get("adaptation", 0) < -0.4:
        lines.append("Spike-frequency adaptation likely reduced (tonic/burst-prone firing).")
    return lines or ["Mixed coupling; measure subthreshold IV, Rin, and f–I curves to resolve sign."]


def _append_msc_handover(
    lines: list[str],
    *,
    spec: PlanningSpec,
    order_table: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    """Two-month MSc project timeline and milestones."""
    lines.extend([
        "",
        "7. MSc STUDENT PROJECT HANDOVER (8 weeks)",
        "-" * 40,
        f"Project: {spec.msc_project_title}",
        "Primary neuromodulator: Acetylcholine (ACh). Secondary: Noradrenaline (NE).",
        "Target cells: V2M L5 IT and L5 ET pyramidal supertypes (11 Allen supertypes total).",
        "",
        "Week 1 — Onboarding & validation",
        "  • Read receptor_excitability.md §1.4 (HCN axis) and §3 (V2M afferent presynaptic logic).",
        "  • Review expression heatmaps and IT vs ET delta figures in this folder.",
        "  • Practice whole-cell in VISpm/VISam/RSPagl slices; note layer and morphology (IT vs ET).",
        "",
        "Weeks 2–4 — ACh phase (PRIMARY)",
        "  • Order: pirenzepine (M1), methoctramine/AF-DX116 (M2), oxotremorine-M, PNU-282987/MLA (α7),",
        "    nicotine or cytisine + DHβE (α4β2), mecamylamine (broad nAChR backup).",
        "  • Protocol: baseline Rin, rheobase, f–I (500 ms steps), subthreshold ZAP (HCN vs M-resonance).",
        "  • Apply low then high ACh (or selective agonists); compare L5 IT vs L5 ET within V2M.",
        "  • Success: M1/M3 signature (↓adaptation, ↑Rin) in supertypes with high Chrm1/Chrm3.",
        "",
        "Weeks 5–6 — ACh antagonist dissection",
        "  • Receptor-blocker experiments on best-responding cells from weeks 2–4.",
        "  • Test whether α7 (MLA) vs α4β2 (DHβE) separates IT vs ET nicotinic contributions.",
        "",
        "Weeks 7–8 — NE phase (SECONDARY)",
        "  • Order: clonidine/dexmedetomidine + atipamezole (α2), phenylephrine/prazosin (α1),",
        "    dobutamine/betaxolol (β1).",
        "  • Same intrinsic protocol; emphasise Rin and I_h (α2 ↓I_h vs β1 ↑I_h opposition).",
        "  • Compare NE low-tone (α2-dominated) vs high-tone (α1+β1) predictions to data.",
        "",
        "Deliverables",
        "  • Cell table: supertype (morphology), subregion, receptor hits expected, pharmacology applied.",
        "  • Summary figure: IT vs ET × ACh/NE effect sizes on Rin, I_h proxy, adaptation.",
        "  • Short methods paragraph linking back to this expression evidence folder.",
        "",
    ])
    if spec.compare_groups:
        ga, gb = spec.compare_groups
        lines.append(f"Key IT vs ET expression differences to test ({gb} − {ga}):")
        for family in spec.family_order:
            sub = summary[summary["family"] == family]
            genes = sorted(sub["gene"].unique())
            for gene in genes:
                d = _compare_group_delta(
                    summary, gene, ga, gb,
                    brain_area=spec.brain_areas[0] if spec.brain_areas else None,
                )
                if np.isfinite(d) and abs(d) > 0.3:
                    lines.append(f"  {FAMILY_LABEL.get(family, family)} {gene}: Δ={d:+.2f} log2 CPM+1")
        lines.append("")


def write_experiment_plan(
    *,
    summary: pd.DataFrame,
    scenario_summary: pd.DataFrame,
    order_table: pd.DataFrame,
    scenario_table: pd.DataFrame,
    synthesis_run: Path,
    output_path: Path | str,
    cell_type_level: str,
    spec: PlanningSpec | None = None,
) -> Path:
    """Write plain-text experiment plan for NE and ACh."""
    spec = spec or DEFAULT_VISp_V2M_SPEC
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fam_lbl = {"noradrenaline": "Noradrenaline (NE)", "acetylcholine": "Acetylcholine (ACh)"}
    lines: list[str] = [
        spec.report_title,
        "=" * 72,
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"Evidence source: {synthesis_run}",
        f"Cell type level: {cell_type_level}",
        spec.region_label,
        "Expression tiers: high (≥2 independent measured datasets), medium (1 measured).",
        COUPLING_LEGEND + ".",
        "",
        "1. EXECUTIVE SUMMARY",
        "-" * 40,
    ]

    if spec.include_msc_handover:
        lines.append(f"Project focus: {spec.msc_project_title}")
        lines.append("Modulator priority: (1) ACh, (2) NE.")
        lines.append("")

    for family in spec.family_order:
        top = order_table[
            order_table["gene"].isin(summary.loc[summary["family"] == family, "gene"].unique())
        ].head(8)
        high_genes = summary[
            (summary["family"] == family) & (summary["confidence_tier"] == "high")
        ]["gene"].unique()
        scope = "V2M L5 IT/ET supertypes" if spec.per_supertype else "L2/3–L5 coarse types"
        lines.append(
            f"{fam_lbl[family]}: {len(high_genes)} genes at high confidence in {scope}."
        )
        if not top.empty:
            def _label(r: pd.Series) -> str:
                ct = r.get("cell_type") or r["top_cell_type"]
                return f"{r.gene} ({ct})"

            order_str = ", ".join(_label(r) for _, r in top.head(5).iterrows())
            lines.append(f"  Top ordering priorities: {order_str}")

    lines.extend([
        "",
        "2. COMPOUNDS TO ORDER (receptor-selective pharmacology)",
        "-" * 40,
        "Priority = high/medium expression × postsynaptic intrinsic relevance.",
        "Order ACh toolkit before NE toolkit (project priority).",
        "",
    ])

    for family in spec.family_order:
        lines.append(f"--- {fam_lbl[family]} ---")
        sub = order_table[
            order_table["gene"].isin(
                summary.loc[summary["family"] == family, "gene"].unique(),
            )
        ].head(18)
        for _, r in sub.iterrows():
            ct = r.get("cell_type") or r["top_cell_type"]
            lines.append(
                f"  [{r['confidence_tier'].upper()}] {r['receptor']} ({r['gene']}) — "
                f"{r['coarse_type']}, {r['brain_area']} | {ct}"
            )
            lines.append(f"      Agonists: {r['agonists']}")
            lines.append(f"      Antagonists: {r['antagonists']}")
        lines.append("")

    lines.extend([
        "3. PER-RECEPTOR INTRINSIC EXCITABILITY PREDICTIONS",
        "-" * 40,
        "'Excitability' axis in figures = heuristic rest-depolarisation score; see rin/ih/mcurrent for α2-like effects.",
        "See receptor_excitability.md for full mechanisms.",
        "",
    ])
    shown: set[str] = set()
    for _, r in order_table.head(30).iterrows():
        key = str(r["gene"])
        if key in shown:
            continue
        shown.add(key)
        meta = get_pharmacology(key)
        if meta is None:
            continue
        lines.append(
            f"  {meta.receptor} ({key}) [{meta.coupling}]: {meta.intrinsic_effect}"
        )
        lines.append(
            f"      Rin {meta.rin_effect} | I_h {meta.ih_effect} | "
            f"M-current {meta.mcurrent_effect} | Firing: {meta.firing_effect}"
        )

    lines.extend([
        "",
        "4. COMBINED TRANSMITTER APPLICATION (NE or ACh bath)",
        "-" * 40,
        "Low vs high tone uses relative receptor affinities × expression weights.",
        "NE low: α2A/α2C dominant (Gi, ↓I_h, ↑Rin). NE high: α1 + β1 engage (Gq + Gs/Gβγ ↑I_h).",
        "ACh low: M2/M4 presynaptic Gi (↓ release). ACh high: M1/M3 ↓M-current + nAChR depolarisation.",
        "",
    ])

    coarse_types = list(spec.coarse_types or ["L2/3 IT", "L5 IT", "L5 ET"])
    brain_areas = list(spec.brain_areas)
    for family in spec.family_order:
        lines.append(f"--- {fam_lbl[family]} scenarios ---")
        for coarse in coarse_types:
            for area in brain_areas:
                for conc in ("low", "high"):
                    scen = compute_transmitter_scenario(
                        scenario_summary, coarse_type=coarse, brain_area=area,
                        family=family, concentration=conc,
                    )
                    if not scen["active_genes"]:
                        continue
                    lines.append(
                        f"  {coarse} · {area} · {conc.upper()} {scen['transmitter']}: "
                        f"genes={', '.join(scen['active_genes'])}"
                    )
                    fr = scen["coupling_fraction"]
                    lines.append(
                        "      Coupling mix: "
                        + ", ".join(f"{k}={v:.0%}" for k, v in fr.items() if v > 0.05)
                    )
                    for bullet in _interpret_axis_scores(scen["axis_scores"]):
                        lines.append(f"      → {bullet}")
        lines.append("")

    lines.extend([
        "5. SUGGESTED EPHYS READOUTS",
        "-" * 40,
        "  • Input resistance (Rin) — α2 ↓I_h vs β1 ↑I_h vs M1 ↓M-current",
        "  • Subthreshold resonance (ZAP or chirp) — HCN vs M-resonance axes",
        "  • f–I curve / first-spike latency — adaptation loss under M1/M3 or α1",
        "  • Rheobase and AP threshold — Gq vs Gi balance",
        "  • Synaptic isolation (CNQX/D-AP5/GABAzine) for pure intrinsic effects",
        "",
        "6. CAVEATS",
        "-" * 40,
        "  • Expression from Allen MERFISH (+ Vizgen/Zhuang concordance); some MERFISH imputed.",
        "  • Supertype labels are transcriptomic — validate IT vs ET by morphology/projection.",
        "  • Presynaptic receptors (M2, α2C, α7 presynaptic) alter synaptic drive.",
        "  • Circuit effects (ACh → L1 interneuron nicotinic disinhibition) not in intrinsic axes.",
        "  • Affinity weights are ordinal — treat low/high transmitter as qualitative scenarios.",
    ])
    for caveat in spec.extra_caveats:
        lines.append(f"  • {caveat}")
    lines.append("")

    if spec.include_msc_handover:
        _append_msc_handover(lines, spec=spec, order_table=order_table, summary=summary)

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def run_experiment_planning(
    config: dict[str, Any],
    *,
    synthesis_run: Path | str,
    output_dir: Path | str,
    cell_type_level: str | None = None,
    spec: PlanningSpec | None = None,
) -> dict[str, Any]:
    """End-to-end: load evidence, build tables, figures, and text plan."""
    spec = spec or DEFAULT_VISp_V2M_SPEC
    level = cell_type_level or config["cell_type_level"]
    synthesis_run = Path(synthesis_run)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    evidence = load_evidence_table(synthesis_run)
    sub = filter_planning_evidence(evidence, config, spec=spec)
    summary = aggregate_expression_summary(sub, per_supertype=spec.per_supertype)
    order_table = build_receptor_order_table(summary)
    scenario_summary = (
        coarse_summary_for_scenarios(summary) if spec.per_supertype else summary
    )
    scenario_table = build_scenario_table(
        scenario_summary,
        coarse_types=list(spec.coarse_types or ["L2/3 IT", "L5 IT", "L5 ET"]),
        brain_areas=list(spec.brain_areas),
        families=spec.family_order,
    )

    paths: dict[str, Path] = {}
    fig_idx = 1
    for family in spec.family_order:
        fam = FAMILY_LABEL[family]

        def _fig(name: str) -> Path:
            nonlocal fig_idx
            path = output_dir / f"fig{fig_idx:02d}_{name}.png"
            fig_idx += 1
            return path

        paths[f"heatmap_{fam}"] = plot_expression_heatmap(
            summary, family, _fig(f"expression_{fam.lower()}"), config=config, spec=spec,
        )
        if spec.per_supertype:
            paths[f"dot_{fam}"] = plot_subtype_dot_comparison(
                summary, family, _fig(f"supertype_dots_{fam.lower()}"), config=config,
            )
        paths[f"priority_{fam}"] = plot_priority_bars(
            order_table, family,
            _fig(f"priority_order_{fam.lower()}"),
            config=config,
            title_suffix=" (V2M L5)" if spec.per_supertype else "",
        )
        if spec.include_region_delta:
            paths[f"delta_{fam}"] = plot_region_delta(
                summary, family, _fig(f"visp_v2m_delta_{fam.lower()}"), config=config,
            )
        elif spec.compare_groups:
            paths[f"delta_{fam}"] = plot_it_et_delta(
                summary, family, _fig(f"it_et_delta_{fam.lower()}"),
                config=config, spec=spec,
            )
        paths[f"coupling_{fam}"] = plot_coupling_scenarios(
            scenario_table, family, _fig(f"coupling_scenario_{fam.lower()}"), config=config,
        )
        paths[f"axis_{fam}"] = plot_axis_scenarios(
            scenario_table, family, _fig(f"intrinsic_axes_{fam.lower()}"), config=config,
        )

    summary.to_csv(output_dir / "expression_summary.csv", index=False)
    order_table.to_csv(output_dir / "pharmacology_order_table.csv", index=False)
    scenario_table.to_csv(output_dir / "transmitter_scenario_table.csv", index=False)
    paths["plan"] = write_experiment_plan(
        summary=summary,
        scenario_summary=scenario_summary,
        order_table=order_table,
        scenario_table=scenario_table,
        synthesis_run=synthesis_run,
        output_path=output_dir / spec.plan_filename,
        cell_type_level=level,
        spec=spec,
    )
    return {
        "summary": summary,
        "order_table": order_table,
        "scenario_table": scenario_table,
        "paths": {k: v for k, v in paths.items() if v is not None},
        "synthesis_run": synthesis_run,
        "output_dir": output_dir,
        "spec": spec,
    }
