"""NE/ACh experiment planning from synthesis evidence tables.

Loads prior synthesis outputs, ranks pharmacology targets, predicts combined
transmitter effects at low vs high concentration, writes figures and a text plan.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import matplotlib.pyplot as plt
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
AXIS_LABELS = ["excitability", "rin", "ih", "mcurrent", "adaptation"]


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
    families: tuple[str, ...] = ("noradrenaline", "acetylcholine"),
    brain_areas: tuple[str, ...] = ("VISp", "V2M"),
    min_tier: str = "medium",
) -> pd.DataFrame:
    """High/medium NE or ACh rows for VISp/V2M pyramidal-filtered cell types."""
    allowed_tiers = {"high", "medium"} if min_tier == "medium" else {"high"}
    sub = evidence[
        evidence["family"].isin(families)
        & evidence["brain_area"].isin(brain_areas)
        & evidence["confidence_tier"].isin(allowed_tiers)
    ].copy()

    name_filter = config.get("cell_type_name_filter") or []
    if name_filter:
        mask = pd.Series(False, index=sub.index)
        for token in name_filter:
            mask |= sub["cell_type"].str.contains(token, na=False, regex=False)
        sub = sub[mask]

    sub["coarse_type"] = sub["cell_type"].map(lambda s: infer_coarse_type(str(s)))
    sub = sub[sub["coarse_type"].notna()].copy()

    expr_vals: list[float] = []
    expr_ds: list[str | None] = []
    for _, row in sub.iterrows():
        val, ds, _src = pick_expression_value(row, config)
        expr_vals.append(val)
        expr_ds.append(ds)
    sub["expression"] = expr_vals
    sub["expression_dataset"] = expr_ds
    return sub


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
    return out.sort_values(
        ["coarse_type", "brain_area", "priority_score"],
        ascending=[True, True, False],
    )


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
) -> Path:
    """Coarse type × gene heatmap with VISp and V2M side-by-side columns."""
    sub = summary[summary["family"] == family].copy()
    if sub.empty:
        raise ValueError(f"No summary rows for family {family!r}")

    coarse_order = ["L2/3 IT", "L5 IT", "L5 ET"]
    genes = sorted(sub["gene"].unique(), key=lambda g: (
        -sub.loc[sub["gene"] == g, "expression"].max(),
        g,
    ))
    areas = ["VISp", "V2M"]
    index = pd.MultiIndex.from_product([coarse_order, areas], names=["coarse_type", "brain_area"])
    mat = pd.DataFrame(index=index, columns=genes, dtype=float)

    for _, row in sub.iterrows():
        key = (row["coarse_type"], row["brain_area"])
        if key not in mat.index:
            continue
        mat.loc[key, row["gene"]] = row["expression"]

    row_labels = [f"{c} · {a}" for c, a in mat.index]
    fig_w = max(10, 0.45 * len(genes) + 4)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
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
    fam_lbl = FAMILY_LABEL.get(family, family)
    ax.set_title(f"{fam_lbl} receptor expression — coarse pyramidal types (VISp vs V2M)")
    ax.set_xlabel("Gene")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=config["output"].get("dpi", 150), bbox_inches="tight")
    plt.close(fig)
    return path


def plot_priority_bars(
    order_table: pd.DataFrame,
    family: str,
    output_path: Path | str,
    *,
    config: dict[str, Any],
    top_n: int = 12,
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
    plot_df["label"] = plot_df.apply(
        lambda r: f"{r['gene']} ({r['coarse_type']}, {r['brain_area']})", axis=1,
    )

    fig, ax = plt.subplots(figsize=(10, max(5, 0.35 * len(plot_df))))
    colors = plot_df["coupling"].map(COUPLING_COLORS).fillna("#999999")
    ax.barh(plot_df["label"], plot_df["priority_score"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("Priority score (tier + expression + coupling)")
    fam_lbl = FAMILY_LABEL.get(family, family)
    ax.set_title(f"{fam_lbl} — pharmacology ordering candidates")
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


def write_experiment_plan(
    *,
    summary: pd.DataFrame,
    order_table: pd.DataFrame,
    scenario_table: pd.DataFrame,
    synthesis_run: Path,
    output_path: Path | str,
    cell_type_level: str,
) -> Path:
    """Write plain-text experiment plan for NE and ACh."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fam_lbl = {"noradrenaline": "Noradrenaline (NE)", "acetylcholine": "Acetylcholine (ACh)"}
    lines: list[str] = [
        "NE & ACh RECEPTOR EXPRESSION — EXPERIMENT PLAN",
        "=" * 72,
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"Evidence source: {synthesis_run}",
        f"Cell type level: {cell_type_level}",
        "Regions: VISp (V1) vs V2M (VISpm + VISam + RSPagl unweighted rollup)",
        "Expression tiers: high (≥2 independent measured datasets), medium (1 measured).",
        "",
        "1. EXECUTIVE SUMMARY",
        "-" * 40,
    ]

    for family in ("noradrenaline", "acetylcholine"):
        fam_short = FAMILY_LABEL[family]
        top = order_table[order_table["gene"].isin(summary.loc[summary["family"] == family, "gene"])].head(8)
        high_genes = summary[
            (summary["family"] == family) & (summary["confidence_tier"] == "high")
        ]["gene"].unique()
        lines.append(
            f"{fam_lbl[family]}: {len(high_genes)} genes at high confidence across "
            f"L2/3–L5 coarse types in VISp/V2M."
        )
        if not top.empty:
            order_str = ", ".join(
                f"{r.gene} ({r.coarse_type}, {r.brain_area})" for _, r in top.head(5).iterrows()
            )
            lines.append(f"  Top ordering priorities: {order_str}")

    lines.extend([
        "",
        "2. COMPOUNDS TO ORDER (receptor-selective pharmacology)",
        "-" * 40,
        "Priority = high/medium expression in target pyramidal class × postsynaptic intrinsic relevance.",
        "",
    ])

    for family in ("noradrenaline", "acetylcholine"):
        lines.append(f"--- {fam_lbl[family]} ---")
        sub = order_table[
            order_table["gene"].isin(
                summary.loc[summary["family"] == family, "gene"].unique(),
            )
        ].head(15)
        for _, r in sub.iterrows():
            lines.append(
                f"  [{r['confidence_tier'].upper()}] {r['receptor']} ({r['gene']}) — "
                f"{r['coarse_type']}, {r['brain_area']} | top supertype: {r['top_cell_type']}"
            )
            lines.append(f"      Agonists: {r['agonists']}")
            lines.append(f"      Antagonists: {r['antagonists']}")
        lines.append("")

    lines.extend([
        "3. PER-RECEPTOR INTRINSIC EXCITABILITY PREDICTIONS",
        "-" * 40,
        "Columns: Rin, I_h, M-current, firing. See receptor_excitability.md for mechanisms.",
        "",
    ])
    shown = set()
    for _, r in order_table.head(25).iterrows():
        key = r["gene"]
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
        "NE low: α2A/α2C dominant (Gi, ↓I_h, ↑Rin). NE high: α1 + β1 engage (Gq depolarisation + Gs/Gβγ ↑I_h).",
        "ACh low: M2/M4 presynaptic Gi (↓ release). ACh high: M1/M3 ↓M-current + nAChR depolarisation.",
        "",
    ])

    coarse_types = ["L2/3 IT", "L5 IT", "L5 ET"]
    for family in ("noradrenaline", "acetylcholine"):
        lines.append(f"--- {fam_lbl[family]} scenarios ---")
        for coarse in coarse_types:
            for area in ("VISp", "V2M"):
                for conc in ("low", "high"):
                    scen = compute_transmitter_scenario(
                        summary, coarse_type=coarse, brain_area=area,
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
        "  • Subthreshold resonance (ZAP or chirp) — HCN vs M-current axes",
        "  • f–I curve / first-spike latency — adaptation loss under M1/M3 or α1",
        "  • Rheobase and AP threshold — Gq vs Gi balance",
        "  • Optional: synaptic isolation (CNQX/D-AP5/GABAzine) for pure intrinsic effects",
        "",
        "6. CAVEATS",
        "-" * 40,
        "  • Expression from Allen MERFISH (+ Vizgen/Zhuang concordance); some MERFISH imputed.",
        "  • Supertype/coarse-type aggregation; validate on your recorded cell class.",
        "  • Presynaptic receptors (M2, α2C, α7 presynaptic) alter synaptic drive — use isolated intrinsic protocols.",
        "  • Circuit effects (e.g. ACh → interneuron nicotinic disinhibition) not captured in intrinsic table.",
        "  • Affinity weights are ordinal, not fitted Kd — treat low/high NE/ACh as qualitative scenarios.",
        "",
    ])

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def run_experiment_planning(
    config: dict[str, Any],
    *,
    synthesis_run: Path | str,
    output_dir: Path | str,
    cell_type_level: str | None = None,
) -> dict[str, Any]:
    """End-to-end: load evidence, build tables, figures, and text plan."""
    level = cell_type_level or config["cell_type_level"]
    synthesis_run = Path(synthesis_run)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    evidence = load_evidence_table(synthesis_run)
    sub = filter_planning_evidence(evidence, config)
    summary = aggregate_coarse_summary(sub)
    order_table = build_receptor_order_table(summary)
    scenario_table = build_scenario_table(
        summary,
        coarse_types=["L2/3 IT", "L5 IT", "L5 ET"],
        brain_areas=["VISp", "V2M"],
    )

    paths: dict[str, Path] = {}
    for family in ("noradrenaline", "acetylcholine"):
        fam = FAMILY_LABEL[family]
        paths[f"heatmap_{fam}"] = plot_expression_heatmap(
            summary, family, output_dir / f"fig01_expression_{fam.lower()}.png", config=config,
        )
        paths[f"priority_{fam}"] = plot_priority_bars(
            order_table, family, output_dir / f"fig02_priority_order_{fam.lower()}.png", config=config,
        )
        paths[f"delta_{fam}"] = plot_region_delta(
            summary, family, output_dir / f"fig03_visp_v2m_delta_{fam.lower()}.png", config=config,
        )
        paths[f"coupling_{fam}"] = plot_coupling_scenarios(
            scenario_table, family, output_dir / f"fig04_coupling_scenario_{fam.lower()}.png", config=config,
        )
        paths[f"axis_{fam}"] = plot_axis_scenarios(
            scenario_table, family, output_dir / f"fig05_intrinsic_axes_{fam.lower()}.png", config=config,
        )

    summary.to_csv(output_dir / "coarse_expression_summary.csv", index=False)
    order_table.to_csv(output_dir / "pharmacology_order_table.csv", index=False)
    scenario_table.to_csv(output_dir / "transmitter_scenario_table.csv", index=False)
    paths["plan"] = write_experiment_plan(
        summary=summary,
        order_table=order_table,
        scenario_table=scenario_table,
        synthesis_run=synthesis_run,
        output_path=output_dir / "EXPERIMENT_PLAN_NE_ACh.txt",
        cell_type_level=level,
    )
    return {
        "summary": summary,
        "order_table": order_table,
        "scenario_table": scenario_table,
        "paths": paths,
        "synthesis_run": synthesis_run,
        "output_dir": output_dir,
    }
