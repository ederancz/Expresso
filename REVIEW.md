# Expresso — Code & Science Review

Scope: full-repository review with a focus on **scientific correctness**, oriented
toward the project goal — using multiple openly available data sources to synthesise
and cross-validate the expression of receptor and excitability gene panels in
specific brain areas and cell types, within the Allen CCF + Allen cell-type
taxonomy framework, to support statements of the form *"cell type C in region R
expresses gene set X, implying excitability profile … and neuromodulatory
influences …"*.

Status: **review complete; all approved fixes implemented.** This document records
the original findings and, for each, the resolution now in the codebase. Items are
tagged **[Fixed]**, **[Documented]** (accepted/explained, no code change), or
**[Open]**.

---

## 0. Overall verdict

The architecture is sound and well-suited to the goal:

- Config-driven (no hard-coded genes/regions in logic).
- Memory-safe backed/partial `h5ad` reads with explicit file-handle closing.
- Per-run timestamped output folders with `run_manifest.json` capturing git
  commit + config snapshot (good reproducibility).
- Graceful handling of genes missing from a given dataset.
- Four datasets wired for cross-validation: Allen scRNA (WMB-10Xv3), Allen MERFISH
  (± imputed), Vizgen receptor map, Zhuang MERFISH.

The correctness issues that could change conclusions have been addressed, the
**excitability arm now runs** off the unified config, and the **synthesis /
cross-validation statement step now exists** (`src/synthesis.py` +
`notebooks/05_synthesis.ipynb`). Remaining items are low-severity polish, tracked
in §5–§6.

---

## 1. Neurobiological logic (as understood)

- **Cell type (e.g. L5 ET):** L5 extratelencephalic / pyramidal-tract neurons —
  Allen subclass `L5 ET CTX Glut`. Distinctive intrinsic excitability (HCN/Ih sag,
  large apical dendritic Ca²⁺ spikes, burst-capable) and rich neuromodulation, so
  a sensible target for an excitability + neuromodulation story. **The relevant
  taxonomy level and exact clusters are deliberately left open** — the pipeline must
  operate at `class`/`subclass`/`supertype`/`cluster` and the biologically relevant
  level/clusters will be chosen later.
- **Region (e.g. VISpm):** a fine CCF parcel. Distinguishing VISpm from VISp
  requires spatial data (Allen/Zhuang MERFISH CCF parcellation, or Vizgen via label
  transfer). Dissociated scRNA only carries the coarse `VIS` dissection ROI.
- **Receptor panel → neuromodulatory influences:** receptor families expressed
  imply which neuromodulators act, and Gs/Gi/Gq coupling implies the sign.
- **Excitability panel → intrinsic profile:** Nav/Kv/Cav/HCN/KCa/Kir + auxiliary
  subunits + Ca²⁺-handling genes → firing/integration phenotype. `excitability_genes.md`
  already maps each gene to mechanism and a PROMINENT/PLAUSIBLE tier.
- **Cross-validation:** the same gene × cell-type × region quantity measured in
  independent resources; concordance across modalities is what makes an expression
  call trustworthy.

Scientific weak points to respect: VISpm specificity must come from spatial
datasets (not scRNA); mean-with-zeros needs a companion fraction-expressing; many
excitability/peptide genes are imputed (not measured) in Allen MERFISH.

---

## 2. Scientific correctness (priority)

### 2.1 [Fixed] Vizgen "CPM" normalized over the gene subset, not library size
Previously `load_vizgen_sample` read only the requested receptor columns, then
divided each cell by the row sum across *those* genes — the denominator was "total
receptor counts", not library size. This inflated values (saved cross-ref showed
Vizgen ~9–12 vs Allen ~1–3), broke comparability, and created artifactual
structure for cells expressing few receptors.

**Resolution:** `load_vizgen_sample` now computes per-cell library totals over the
**full real-gene panel** (all non-`Blank*` columns), and
`_log2_cpm_plus_one_with_totals` applies `log2(CPM+1)` for the requested subset
using those totals. No volume division.

### 2.2 [Fixed] Cross-dataset absolute comparison mixes CPM denominators
Panels differ in size (Allen ≈500, Zhuang ≈1122, Vizgen ≈483), so log2(CPM+1)
magnitudes are not directly comparable across datasets. The scatter previously led
with **Pearson r** and a y=x identity line, over-promising comparability.

**Resolution:** `plot_crossref_scatter` now leads with **Spearman ρ** (rank
concordance), reports a measured-only ρ as well, and demotes Pearson r with a
caveat about differing CPM denominators (`_safe_corr` helper).

### 2.3 [Fixed] kNN label transfer in unnormalized, unscaled space
`transfer_allen_merfish_labels` did Euclidean kNN with Allen log2(CPM+1) as
reference and Vizgen as query, with no per-gene scaling — making the transferred
cell-type/region labels (which underpin any Vizgen-based claim) the least reliable
link.

**Resolution:** features are now z-scored per gene on the shared genes
(`_standardize_fit` fits on the reference, applied to both), and
`_knn_majority_labels_with_confidence` records a neighbour-vote confidence for both
cell type and brain area. `load_vizgen_aggregated` can drop low-confidence cells
via `vizgen_label_transfer_min_confidence`.

### 2.4 [Fixed] Guard against circular cross-validation (imputed vs measured)
Most genes are imputed in Allen MERFISH (saved run: 75/109 receptors; worse for
excitability ion channels). Because Allen imputation uses WMB-10x scRNA,
Allen-imputed ↔ Allen-scRNA agreement is partly circular, whereas Allen-imputed ↔
Zhuang-measured and ↔ Vizgen-measured are independent.

**Resolution:** provenance is tracked end-to-end. Imputed genes are marked in plots
(`*`, hollow markers); the cross-ref reports measured-only metrics separately; and
the synthesis evidence table counts `n_independent_measured_detections` and assigns
confidence tiers that lean on genes measured in ≥2 independent datasets (§4).

### 2.5 [Fixed] Mean-with-zeros is the only summary statistic
`aggregate_scrna_expression` reported only mean log2(CPM+1) over all cells. A
defensible "expresses" call also needs **fraction of cells expressing**.

**Resolution:** aggregation now returns `frac_expressing` and `n_cells` alongside
`mean_expression`; replicate averaging (`aggregate_zhuang_replicates_mean`, also
used for Vizgen) averages `frac_expressing` and sums `n_cells`. The synthesis dot
plot uses size = fraction expressing, colour = mean.

### 2.6 [Documented] scRNA region resolution vs a fine-region claim
Handled correctly in code (pooling warning): notebook 01 cannot support a
VISpm-specific claim (it pools all VIS dissection). Fine-region specificity comes
from Allen MERFISH / Zhuang (native CCF) and Vizgen (label transfer). The synthesis
encodes this split: scRNA → cell-type-level expression; spatial datasets → region
localization (`region_resolved` flag per dataset).

---

## 3. Excitability arm (core; now runnable)

### 3.1 [Fixed] Config key mismatch
`config.py` previously hard-required a top-level `receptors` key, and notebooks
hard-coded `receptor_query_config.yaml`, so the excitability panel could not run.

**Resolution:** a generic `gene_panel` key was introduced (`_parse_gene_panel`),
accepting either a flat `family → genes` map or a nested
`category → family → genes` map; legacy top-level `receptors`/`excitability` keys
still load. Both panels now live in the unified `query_config.yaml` under
`gene_panel` (categories `receptors` + `excitability`), and all notebooks point at
it. The old per-panel YAMLs (`receptor_query_config.yaml`,
`excitability_query_config.yaml`) have been removed.

### 3.2 [Open] Excitability genes mostly outside MERFISH panels
Allen MERFISH imputes/lacks most ion-channel genes and Vizgen (a receptor panel)
covers almost none, so for excitability the credible datasets are Allen scRNA (all
genes, cell-type resolution) + Zhuang (region localization, partial coverage), with
Allen-imputed as supporting — exactly what the synthesis confidence tiers capture.
*Optional future enhancement:* lift the PROMINENT/PLAUSIBLE tier from
`excitability_genes.md` into the YAML so it can drive prioritization.

---

## 4. Synthesis + cross-validation — [Fixed / implemented]

Implemented in `src/synthesis.py` and driven by `notebooks/05_synthesis.ipynb`.
Level-agnostic by design: operates at whatever `cell_type_level` is configured
(default `supertype`); "L5 ET × VISpm" is only an illustrative target.

- **Dataset registry & discovery:** `DATASET_SPECS` describes each dataset
  (parquet names, region resolution, imputed support, independence);
  `discover_dataset_parquets` / `gather_dataset_aggregates` locate and load the
  newest per-dataset aggregates from prior runs.
- **Per-gene, per-dataset evidence table** (`build_evidence_table`): one row per
  (cell type × region × gene); per-dataset `mean_expression`, `frac_expressing`,
  `n_cells`, and `source` (measured/imputed for Allen). scRNA is broadcast across
  regions and flagged `region_resolved=False`. A `category` column carries
  receptor/excitability.
- **Cross-validation metrics:** configurable detection flag (`frac ≥ min_frac` AND
  `mean ≥ min_mean`), `n_independent_measured_detections` (+ a separate
  `supporting_imputed_detection`), and a `confidence_tier` (high = measured in ≥2
  independent datasets; medium = 1 measured + Allen-imputed; low = single/imputed).
- **Targeted outputs:** `target_evidence` / `summarize_target` filter to a
  cell type (± region); `statement_scaffold` emits a human-readable claim split by
  category (receptors vs excitability); `plot_evidence_dotplot` renders the dot plot
  (mean × fraction across datasets). Tidy evidence is written to parquet/CSV.

---

## 5. Other correctness & mechanics

- [Fixed] Imputed suffix hardcoding: `_imputed_gene_symbols` now takes `config` and
  uses `get_expression_suffix`, consistent with the imputed data slice.
- [Fixed] Duplicated import: `filter_cell_types_by_name` is no longer re-imported
  inside `combined_heatmap_matrix`.
- [Open] `top_variable_cell_types` is retained as a library helper but no longer
  used by `combined_heatmap_matrix`.
- [Open] `plot_spatial` / `plot_family_spatial_panel` remain unused by notebooks
  (NB02 produces heatmaps) — kept as library extras for the spatial-map story.
- [Open] Cross-ref joins assume identical taxonomy strings across datasets (true for
  Allen/Zhuang/transferred Vizgen); an overlap-count log/assertion would guard
  against silent empty joins.

---

## 6. Documentation & reproducibility hygiene

- [Fixed] Removed stale `cursor_handover.md` links from `README.md` /
  `REPOSITORY_GUIDE.md`.
- [Fixed] `README.md` and `REPOSITORY_GUIDE.md` updated for the unified
  `query_config.yaml`, the nested `gene_panel` schema, corrected notebook statuses
  (NB02 = heatmaps; NB03/NB04 implemented; NB05 added), and the removed per-panel
  YAMLs.
- [Documented] `abc_atlas_access` is pinned to git HEAD; the per-run manifest
  snapshot mitigates this for reproducibility.

---

## 7. Fix plan — completed

Decisions (as approved):
- Vizgen normalization: **full real-gene panel CPM, exclude `Blank*`, no volume**.
- Config: **generic `gene_panel` key**, unified into `query_config.yaml`.
- Cell-type level: **level-agnostic**; biologically relevant level/clusters TBD later.
- Region guardrail on cross-run parquet discovery: **warn** on mismatch.

All items below are implemented:
1. ✅ §2.1 Vizgen full-library normalization.
2. ✅ §3.1 generic `gene_panel` config support (excitability arm unlocked).
3. ✅ §2.5 fraction-expressing + `n_cells` in aggregation.
4. ✅ §2.2 / §2.4 Spearman-first stats + measured-vs-imputed split in cross-ref.
5. ✅ §2.3 z-scored kNN label transfer + confidence reporting.
6. ✅ §4 level-agnostic synthesis module + notebook.
7. ✅ §5 / §6 mechanics cleanup + doc sync (incl. region-mismatch warning guardrail).
