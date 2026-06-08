# Expresso — Code & Science Review

Scope: full-repository review with a focus on **scientific correctness**, oriented
toward the project goal — using multiple openly available data sources to synthesise
and cross-validate the expression of receptor and excitability gene panels in
specific brain areas and cell types, within the Allen CCF + Allen cell-type
taxonomy framework, to support statements of the form *"cell type C in region R
expresses gene set X, implying excitability profile … and neuromodulatory
influences …"*.

Status: **review only.** No source code changed by this document. Fix decisions
captured at the end.

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

However, for the **specific** scientific claim (cell type × region expresses gene
set, cross-validated), there are correctness issues that can change conclusions,
the **excitability arm is currently non-runnable**, and the **final synthesis /
cross-validation statement step does not yet exist**.

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

### 2.1 [Critical] Vizgen "CPM" normalized over the gene subset, not library size
`load_vizgen_sample` reads only the requested receptor columns, then
`_normalize_log2_cpm_plus_one` divides each cell by the row sum across *those*
genes — i.e. the denominator is "total receptor counts", not library size.

Consequences:
- Inflated values (saved cross-ref shows Vizgen ~9–12 vs Allen ~1–3).
- Not CPM; not comparable to Allen/Zhuang.
- Cells expressing few receptors get each one scaled up → artifactual structure.

This corrupts the Allen↔Vizgen cross-validation. **Decision: normalize over the
full real-gene panel (exclude `Blank*` control probes); no volume division.**

### 2.2 [Major] Cross-dataset absolute comparison mixes CPM denominators
Even after 2.1, panels differ (Allen ≈500, Zhuang ≈1122, Vizgen ≈483), so
log2(CPM+1) magnitudes are not directly comparable across datasets. The scatter
currently leads with **Pearson r** and an identity (y=x) line, over-promising
comparability.

Fix: make **Spearman ρ** (rank concordance) the headline; keep Pearson secondary;
drop/annotate the y=x line or compare per-dataset z-scored values.

### 2.3 [Major] kNN label transfer in unnormalized, unscaled space
`transfer_allen_merfish_labels` does Euclidean kNN with Allen log2(CPM+1) as
reference and Vizgen (subset-normalized) as query, no per-gene scaling/PCA. The
resulting cell-type/region labels on Vizgen cells — the cells underpinning a
Vizgen-based claim — are the least reliable link.

Fix: z-score per gene on the shared genes (fit on reference, apply to query)
before kNN; report transfer confidence (neighbour-vote fraction) so low-confidence
cells can be excluded.

### 2.4 [Major, conceptual] Guard against circular cross-validation (imputed vs measured)
Most genes are imputed in Allen MERFISH (saved run: 75/109 receptors; worse for
excitability ion channels). Allen imputation uses WMB-10x scRNA, so:
- Allen-imputed ↔ Allen-scRNA agreement is partly circular (not independent).
- Allen-imputed ↔ Zhuang-measured and ↔ Vizgen-measured are independent.

The code marks imputed genes (`*`, hollow markers) — good — but the
cross-validation **metric** should be reported **separately for measured vs
imputed** Allen genes, and the final claim should lean on genes measured in ≥2
independent datasets.

### 2.5 [Moderate] Mean-with-zeros is the only summary statistic
`aggregate_scrna_expression` reports mean log2(CPM+1) over all cells (zeros
included). A defensible "expresses" call also needs **fraction of cells
expressing** (detection rate). Add fraction-expressing to the aggregation; the
standard artifact is a dot plot (color = mean-in-expressing, size = fraction).

### 2.6 [Moderate] scRNA region resolution vs a fine-region claim
Handled correctly in code (pooling warning), but to state explicitly: notebook 01
cannot support a VISpm-specific claim (pools all VIS dissection). Fine-region
specificity must come from Allen MERFISH / Zhuang (native CCF) and Vizgen (label
transfer). The synthesis encodes this: scRNA → cell-type-level expression; spatial
datasets → region localization.

---

## 3. Excitability arm (core; currently non-runnable)

### 3.1 [Blocker] Config key mismatch
`config.py` hard-requires top-level `receptors` and builds the gene map from it;
`excitability_query_config.yaml` uses top-level `excitability`, and notebooks
hard-code `receptor_query_config.yaml`. The excitability analysis cannot run.

**Decision: introduce a generic `gene_panel` key and migrate both YAMLs to it.**
All derived keys (`_genes_flat`/`_all_genes`/`_families`) and downstream logic stay
the same; notebooks get a `CONFIG_PATH` switch.

### 3.2 [Moderate] Excitability genes mostly outside MERFISH panels
Expect Allen MERFISH to impute/lack most ion-channel genes; Vizgen (a receptor
panel) covers almost none. For excitability the credible datasets are Allen scRNA
(all genes, cell-type resolution) + Zhuang (region localization, partial coverage),
with Allen-imputed as supporting. Consider lifting the PROMINENT/PLAUSIBLE tier
from `excitability_genes.md` comments into the YAML so it can drive prioritization.

---

## 4. Design spec — synthesis + cross-validation (to build after approval)

Level-agnostic by design: operate at whatever `cell_type_level` is configured;
no hardcoded cell type. "L5 ET × VISpm" is only an illustrative target.

- **Target selection (optional config block):** `target: {cell_type: <name or null>,
  region: <acronym or null>}`. When null, produce the full table; when set, produce
  the focused claim. Selection happens at the configured level.
- **Per-gene, per-dataset evidence table:** one row per gene; per-dataset columns
  for `mean_expr`, `frac_expressing` (new, §2.5), `n_cells`, `source`
  (measured/imputed for Allen), and per-dataset z-score / within-level rank
  (normalization-robust, §2.2).
- **Cross-validation metrics:**
  - Detection concordance: gene "expressed" in dataset D if `frac_expressing ≥ f`
    AND `mean ≥ m`; count independent datasets detecting it (Allen-imputed counted
    separately, §2.4).
  - Rank concordance: Spearman ρ across cell types between dataset pairs within the
    target region; measured-only and measured+imputed variants.
  - Per-gene confidence tier: high = measured-detected in ≥2 independent datasets;
    medium = 1 measured + Allen-imputed; low = imputed-only / single dataset.
- **Outputs:** tidy `*_evidence.parquet` + CSV; focused dot plot (genes × dataset);
  concordance summary figure; an auto-generated statement scaffold combining a
  curated receptor-family → neuromodulator/sign map and channel-family →
  biophysical-role map (the latter already in `excitability_genes.md`).

---

## 5. Other correctness & mechanics (medium/low)

- [Low] Imputed suffix hardcoded: `_imputed_gene_symbols` reads `imputed/log2`
  regardless of `expression_unit`, while the imputed data slice uses the configured
  suffix. Inconsistent for `expression_unit: raw`.
- [Low] Dead/duplicated helpers: `top_variable_cell_types` no longer used by
  `combined_heatmap_matrix`; `filter_cell_types_by_name` imported at module top and
  inside functions.
- [Low] `plot_spatial` / `plot_family_spatial_panel` are unused by notebooks
  (NB02 produces heatmaps). Either wire spatial maps in (useful for the region
  story) or document them as library extras.
- [Low] Cross-ref join assumes identical taxonomy strings across datasets (true for
  Allen/Zhuang/transferred Vizgen). Add an assertion/log of overlap counts to avoid
  silent empty joins.

---

## 6. Documentation & reproducibility hygiene (low)

- `cursor_handover.md` is deleted but still linked from `README.md` and
  `REPOSITORY_GUIDE.md` — remove links or restore.
- `REPOSITORY_GUIDE.md` is stale: ~50-gene/8-family panel and an "adrenergic"
  family (config uses `noradrenaline`); NB02 described as "spatial scatter maps"
  (it's heatmaps); NB04 described as a "stub" (fully implemented).
- `abc_atlas_access` pinned to git HEAD; manifest snapshot mitigates this for
  reproducibility.

---

## 7. Fix plan (decisions captured)

Decisions:
- Vizgen normalization: **full real-gene panel CPM, exclude `Blank*`, no volume**.
- Config: **generic `gene_panel` key**, migrate both YAMLs.
- Cell-type level: **level-agnostic**; biologically relevant level/clusters TBD later.

Proposed order:
1. §2.1 Vizgen full-library normalization (cross-val correctness blocker).
2. §3.1 generic `gene_panel` config support (unlocks excitability arm).
3. §2.5 add fraction-expressing to aggregation.
4. §2.2 / §2.4 Spearman-first stats + measured-vs-imputed split in cross-ref.
5. §2.3 z-scored kNN label transfer + confidence reporting.
6. §4 build the level-agnostic synthesis module + notebook.
7. §5 / §6 cleanup + doc sync.
