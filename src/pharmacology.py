"""Curated pharmacology and intrinsic-excitability metadata for NE and ACh receptors.

Used by ``experiment_planning`` to link expression evidence to ephys predictions
and compound ordering. Affinity weights are relative (0–1), not absolute Kd values.
See ``receptor_excitability.md`` for mechanism narrative and references.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

Coupling = Literal["Gi", "Gq", "Gs", "ionotropic"]
Concentration = Literal["low", "high"]
Transmitter = Literal["NE", "ACh"]

TIER_ORDER = {"high": 0, "medium": 1, "low": 2, "none": 3}


@dataclass(frozen=True)
class ReceptorPharmacology:
    gene: str
    receptor: str
    family: str
    coupling: Coupling
    intrinsic_effect: str
    rin_effect: str
    ih_effect: str
    mcurrent_effect: str
    firing_effect: str
    agonists: tuple[str, ...]
    antagonists: tuple[str, ...]
    notes: str
    ne_low_weight: float = 0.0
    ne_high_weight: float = 0.0
    ach_low_weight: float = 0.0
    ach_high_weight: float = 0.0


# Relative activation weights at low vs high transmitter (literature order-of-magnitude).
# NE: α2 high-affinity → dominant at tonic/low LC; α1/β1 engage at burst/high LC.
# ACh: M2/M4 high-affinity presynaptic; M1/M3 and nAChR at elevated cholinergic tone.
RECEPTOR_PHARMACOLOGY: dict[str, ReceptorPharmacology] = {
    "Adra1a": ReceptorPharmacology(
        gene="Adra1a", receptor="α1A-AR", family="noradrenaline", coupling="Gq",
        intrinsic_effect="Depolarisation; ↓M-current/SK; facilitates dendritic Ca²+/NMDA spikes",
        rin_effect="↑", ih_effect="—", mcurrent_effect="↓", firing_effect="↑ excitability, ↓ adaptation",
        agonists=("Phenylephrine (α1)", "A61603 (α1A selective)"),
        antagonists=("Prazosin (α1)", "Tamsulosin (α1A)"),
        notes="Postsynaptic Gq; also presynaptic facilitation of glutamate release on afferents.",
        ne_low_weight=0.15, ne_high_weight=0.85,
    ),
    "Adra1b": ReceptorPharmacology(
        gene="Adra1b", receptor="α1B-AR", family="noradrenaline", coupling="Gq",
        intrinsic_effect="Depolarisation; ↓M-current/SK; persistent firing with α2 co-activation",
        rin_effect="↑", ih_effect="—", mcurrent_effect="↓", firing_effect="↑ excitability",
        agonists=("Phenylephrine (α1)",),
        antagonists=("Prazosin (α1)",),
        notes="Often highest α1 transcript in cortex; strong Gq axis.",
        ne_low_weight=0.20, ne_high_weight=1.0,
    ),
    "Adra1d": ReceptorPharmacology(
        gene="Adra1d", receptor="α1D-AR", family="noradrenaline", coupling="Gq",
        intrinsic_effect="Similar α1 Gq depolarisation; lower cortical abundance",
        rin_effect="↑", ih_effect="—", mcurrent_effect="↓", firing_effect="↑ excitability",
        agonists=("Phenylephrine (α1)",),
        antagonists=("Prazosin (α1)", "BMY7378 (α1D)"),
        notes="Variable across datasets; include if medium+ tier in target cell.",
        ne_low_weight=0.10, ne_high_weight=0.75,
    ),
    "Adra2a": ReceptorPharmacology(
        gene="Adra2a", receptor="α2A-AR", family="noradrenaline", coupling="Gi",
        intrinsic_effect="↓I_h (PLC-PKC) → ↑Rin; enhanced temporal summation of distal EPSPs",
        rin_effect="↑↑", ih_effect="↓", mcurrent_effect="—", firing_effect="↑ burst responsiveness",
        agonists=("Clonidine", "Dexmedetomidine (α2 selective)"),
        antagonists=("Yohimbine", "Atipamezole", "Idazoxan"),
        notes="LC autoreceptor; postsynaptic α2A on deep L3/L5 pyramidal cells (PFC data).",
        ne_low_weight=1.0, ne_high_weight=0.75,
    ),
    "Adra2b": ReceptorPharmacology(
        gene="Adra2b", receptor="α2B-AR", family="noradrenaline", coupling="Gi",
        intrinsic_effect="Predicted Gi/GIRK; peripheral-dominant in adults",
        rin_effect="↑", ih_effect="↓", mcurrent_effect="—", firing_effect="uncertain in cortex",
        agonists=("Clonidine (partial α2)",),
        antagonists=("Yohimbine",),
        notes="Low cortical pyramidal relevance unless expressed.",
        ne_low_weight=0.25, ne_high_weight=0.35,
    ),
    "Adra2c": ReceptorPharmacology(
        gene="Adra2c", receptor="α2C-AR", family="noradrenaline", coupling="Gi",
        intrinsic_effect="NE terminal autoreceptor at low [NE]; heteroreceptor ↓ glutamate release",
        rin_effect="↑", ih_effect="↓", mcurrent_effect="—", firing_effect="presynaptic-dominant",
        agonists=("Clonidine", "Dexmedetomidine"),
        antagonists=("Yohimbine", "Atipamezole"),
        notes="Often medium tier; presynaptic gate on NE release and afferents.",
        ne_low_weight=0.90, ne_high_weight=0.70,
    ),
    "Adrb1": ReceptorPharmacology(
        gene="Adrb1", receptor="β1-AR", family="noradrenaline", coupling="Gs",
        intrinsic_effect="↑I_h via Gβγ → depolarisation, ↓Rin, faster membrane dynamics",
        rin_effect="↓", ih_effect="↑", mcurrent_effect="—", firing_effect="↑ subthreshold resonance",
        agonists=("Dobutamine (β1)", "Isoproterenol (β non-selective)"),
        antagonists=("Betaxolol (β1)", "Propranolol (β non-selective)"),
        notes="Key HCN-axis receptor; often highest NE confidence gene in VISp/V2M L2/3–L5.",
        ne_low_weight=0.05, ne_high_weight=1.0,
    ),
    "Adrb2": ReceptorPharmacology(
        gene="Adrb2", receptor="β2-AR", family="noradrenaline", coupling="Gs",
        intrinsic_effect="Predicted ↑cAMP → ↑I_h; lower cortical pyramidal than β1",
        rin_effect="↓", ih_effect="↑", mcurrent_effect="—", firing_effect="↑ excitability (predicted)",
        agonists=("Isoproterenol", "Salbutamol (β2)"),
        antagonists=("Propranolol", "ICI-118551 (β2)"),
        notes="Order if medium+ tier; often imputed-only in MERFISH.",
        ne_low_weight=0.05, ne_high_weight=0.55,
    ),
    "Adrb3": ReceptorPharmacology(
        gene="Adrb3", receptor="β3-AR", family="noradrenaline", coupling="Gs",
        intrinsic_effect="Negligible CNS pyramidal role",
        rin_effect="—", ih_effect="—", mcurrent_effect="—", firing_effect="—",
        agonists=("—",),
        antagonists=("—",),
        notes="Peripheral-dominant; deprioritise unless high tier.",
        ne_low_weight=0.0, ne_high_weight=0.05,
    ),
    "Chrm1": ReceptorPharmacology(
        gene="Chrm1", receptor="M1 mAChR", family="acetylcholine", coupling="Gq",
        intrinsic_effect="↓M-current, ↓SK, ↓Kir2; depolarisation; loss of spike-frequency adaptation",
        rin_effect="↑", ih_effect="—", mcurrent_effect="↓↓", firing_effect="↑ excitability, bursting",
        agonists=("Oxotremorine-M (M1/M3)", "Xanomeline (M1/M4 biased)"),
        antagonists=("Pirenzepine (M1)", "Telenzepine (M1)"),
        notes="Dominant cortical muscarinic on pyramidal cells; top ACh hit in VISp/V2M.",
        ach_low_weight=0.25, ach_high_weight=1.0,
    ),
    "Chrm2": ReceptorPharmacology(
        gene="Chrm2", receptor="M2 mAChR", family="acetylcholine", coupling="Gi",
        intrinsic_effect="Weak postsynaptic GIRK; shapes biphasic ACh response",
        rin_effect="↓", ih_effect="↓", mcurrent_effect="—", firing_effect="presynaptic-dominant",
        agonists=("Oxotremorine-M", "Methacholine"),
        antagonists=("AF-DX 116 (M2)", "Methoctramine (M2)"),
        notes="High-affinity; suppresses glutamate release at IC-type cortical terminals.",
        ach_low_weight=1.0, ach_high_weight=0.55,
    ),
    "Chrm3": ReceptorPharmacology(
        gene="Chrm3", receptor="M3 mAChR", family="acetylcholine", coupling="Gq",
        intrinsic_effect="Gq similar to M1; ↓M-current; dendritic facilitation",
        rin_effect="↑", ih_effect="—", mcurrent_effect="↓", firing_effect="↑ excitability",
        agonists=("Oxotremorine-M", "Bethanechol"),
        antagonists=("Darifenacin (M3)", "4-DAMP (M3)"),
        notes="Often co-expressed with Chrm1 in L2/3 and L5 IT.",
        ach_low_weight=0.20, ach_high_weight=0.90,
    ),
    "Chrm4": ReceptorPharmacology(
        gene="Chrm4", receptor="M4 mAChR", family="acetylcholine", coupling="Gi",
        intrinsic_effect="Weak postsynaptic GIRK; presynaptic heteroreceptor on afferents",
        rin_effect="↓", ih_effect="—", mcurrent_effect="—", firing_effect="↓ synaptic drive",
        agonists=("Oxotremorine-M",),
        antagonists=("PD 102807 (M4)",),
        notes="Interneuron-enriched; lower pyramidal postsynaptic than M2.",
        ach_low_weight=0.75, ach_high_weight=0.40,
    ),
    "Chrm5": ReceptorPharmacology(
        gene="Chrm5", receptor="M5 mAChR", family="acetylcholine", coupling="Gq",
        intrinsic_effect="Lowest cortical mAChR; Gq predicted",
        rin_effect="↑", ih_effect="—", mcurrent_effect="↓", firing_effect="↑ (predicted)",
        agonists=("Oxotremorine-M",),
        antagonists=("—",),
        notes="Rarely high tier; low priority for ordering.",
        ach_low_weight=0.05, ach_high_weight=0.20,
    ),
    "Chrna4": ReceptorPharmacology(
        gene="Chrna4", receptor="nAChR α4 (α4β2)", family="acetylcholine", coupling="ionotropic",
        intrinsic_effect="Sustained inward cation current; depolarisation L2–L6 pyramidal",
        rin_effect="↓", ih_effect="—", mcurrent_effect="—", firing_effect="↑ firing (direct + circuit)",
        agonists=("Nicotine", "Cytisine (α4β2 partial)", "Varenicline"),
        antagonists=("Mecamylamine (broad nAChR)", "Dihydro-β-erythroidine (α4β2)"),
        notes="Pair with Chrnb2 for α4β2; circuit disinhibition via L1 interneurons may dominate.",
        ach_low_weight=0.35, ach_high_weight=0.85,
    ),
    "Chrna7": ReceptorPharmacology(
        gene="Chrna7", receptor="nAChR α7", family="acetylcholine", coupling="ionotropic",
        intrinsic_effect="Fast desensitising Ca²+-permeable current; presynaptic glutamate facilitation",
        rin_effect="↓", ih_effect="—", mcurrent_effect="—", firing_effect="phasic ↑ (desensitises)",
        agonists=("PNU-282987 (α7)", "GTS-21"),
        antagonists=("MLA (α7)", "Methyllycaconitine"),
        notes="High affinity but rapid desensitisation; TC afferents nicotinic-facilitated (V2M context).",
        ach_low_weight=0.65, ach_high_weight=0.45,
    ),
    "Chrna2": ReceptorPharmacology(
        gene="Chrna2", receptor="nAChR α2", family="acetylcholine", coupling="ionotropic",
        intrinsic_effect="Marks L5 IT corticocortical subtype; modest direct current",
        rin_effect="↓", ih_effect="—", mcurrent_effect="—", firing_effect="subtype marker",
        agonists=("Nicotine",),
        antagonists=("Mecamylamine",),
        notes="Useful as cell-type marker; lower priority for pharmacology unless high tier.",
        ach_low_weight=0.20, ach_high_weight=0.35,
    ),
    "Chrnb2": ReceptorPharmacology(
        gene="Chrnb2", receptor="nAChR β2 (α4β2)", family="acetylcholine", coupling="ionotropic",
        intrinsic_effect="Obligate β subunit for α4β2; sustained depolarisation",
        rin_effect="↓", ih_effect="—", mcurrent_effect="—", firing_effect="↑ excitability",
        agonists=("Nicotine", "Cytisine"),
        antagonists=("Dihydro-β-erythroidine",),
        notes="Co-require Chrna4 expression for α4β2 pharmacology.",
        ach_low_weight=0.30, ach_high_weight=0.80,
    ),
}

for _gene, _extra in {
    "Chrna3": ("nAChR α3", "Autonomic/ganglionic; negligible cortical pyramidal"),
    "Chrna5": ("nAChR α5", "Modulates α4β2; low cortical pyramidal"),
    "Chrna6": ("nAChR α6", "Catecholaminergic-neuron enriched; low pyramidal"),
    "Chrnb3": ("nAChR β3", "Often medium tier in L5 ET; partner subunit context-dependent"),
    "Chrna9": ("nAChR α9", "Cochlear/peripheral; not cortical pyramidal"),
    "Chrna10": ("nAChR α10", "Minimal CNS pyramidal role"),
}.items():
    if _gene not in RECEPTOR_PHARMACOLOGY:
        RECEPTOR_PHARMACOLOGY[_gene] = ReceptorPharmacology(
            gene=_gene, receptor=_extra[0], family="acetylcholine", coupling="ionotropic",
            intrinsic_effect=_extra[1],
            rin_effect="—", ih_effect="—", mcurrent_effect="—", firing_effect="—",
            agonists=("Nicotine (non-selective)",),
            antagonists=("Mecamylamine",),
            notes=_extra[1],
            ach_low_weight=0.15, ach_high_weight=0.25,
        )


def get_pharmacology(gene: str) -> ReceptorPharmacology | None:
    return RECEPTOR_PHARMACOLOGY.get(gene)


def transmitter_weight(gene: str, transmitter: Transmitter, concentration: Concentration) -> float:
    meta = RECEPTOR_PHARMACOLOGY.get(gene)
    if meta is None:
        return 0.0
    if transmitter == "NE":
        return meta.ne_low_weight if concentration == "low" else meta.ne_high_weight
    return meta.ach_low_weight if concentration == "low" else meta.ach_high_weight


def coupling_sign(coupling: Coupling) -> dict[str, float]:
    """Directional contribution to intrinsic excitability axes (−1 inhibitory, +1 excitatory)."""
    if coupling == "Gi":
        return {"excitability": -0.4, "rin": 0.3, "ih": -0.5, "mcurrent": 0.0, "adaptation": 0.2}
    if coupling == "Gq":
        return {"excitability": 0.9, "rin": 0.4, "ih": 0.0, "mcurrent": -0.9, "adaptation": -0.8}
    if coupling == "Gs":
        return {"excitability": 0.5, "rin": -0.4, "ih": 0.8, "mcurrent": 0.0, "adaptation": 0.1}
    return {"excitability": 0.7, "rin": -0.3, "ih": 0.0, "mcurrent": 0.0, "adaptation": -0.2}
