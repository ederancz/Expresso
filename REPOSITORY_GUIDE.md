# Expresso — Repository Deep Dive

A technical guide to the **Expresso** codebase: purpose, architecture, data flows, module APIs, configuration, and operational notes. For quick setup, see [README.md](README.md). For the code/science review and the rationale behind correctness fixes, see [REVIEW.md](REVIEW.md).

---

## Table of contents

1. [What this project does](#what-this-project-does)
2. [Repository layout](#repository-layout)
3. [Architecture overview](#architecture-overview)
4. [Configuration reference](#configuration-reference)
5. [Data sources](#data-sources)
6. [Module reference](#module-reference)
7. [Milestone 1 — scRNA-seq heatmaps](#milestone-1--scrna-seq-heatmaps)
8. [Brain area mapping (scRNA caveat)](#brain-area-mapping-scrna-caveat)
9. [Milestone 2 — MERFISH heatmaps](#milestone-2--merfish-heatmaps)
10. [Milestones 3–4 — Cross-reference](#milestones-34--cross-reference)
11. [Milestone 5 — Synthesis](#milestone-5--synthesis)
12. [Plotting layer](#plotting-layer)
13. [Notebooks](#notebooks)
14. [Outputs and caching](#outputs-and-caching)
15. [Environment and dependencies](#environment-and-dependencies)
16. [Smoke testing](#smoke-testing)
17. [Known limitations and gotchas](#known-limitations-and-gotchas)
18. [Suggested workflows](#suggested-workflows)
19. [Document map](#document-map)

---

## What this project does

Expresso queries **receptor** and **excitability** gene expression across **mouse brain areas** and **cell types** using the [Allen Brain Cell Atlas (ABC Atlas)](https://alleninstitute.github.io/abc_atlas_access/) and complementary MERFISH datasets (Allen, Vizgen, Zhuang). The common reference framework is the **Allen CCF** (brain regions) and the **Allen cell-type taxonomy** (class → subclass → supertype → cluster).

It is designed as a **config-driven analysis pipeline**:

- **No hard-coded genes or regions** in notebook logic — everything comes from [`query_config.yaml`](query_config.yaml).
- **Memory-efficient I/O** — expression matrices are read in backed mode and sliced to only the genes and cells needed.
- **Five implemented milestones:**
  - **M1:** Allen scRNA-seq heatmaps (cell type × brain area)
  - **M2:** Allen MERFISH heatmaps (cell type × CCF brain area; native sub-region resolution)
  - **M3:** Vizgen MERFISH cross-reference vs Allen (label transfer + concordance plots)
  - **M4:** Zhuang MERFISH cross-reference vs Allen (replicate-averaged concordance)
  - **M5:** Cross-dataset synthesis — evidence table, confidence tiers, statement scaffold

The end goal is to support claims of the form: *"cell type C in region R expresses gene set X (receptors + excitability channels), cross-validated across independent datasets, implying defined neuromodulatory influences and intrinsic excitability profiles."*

---

## Repository layout

```
Expresso/
├── README.md                      # Quick start, setup, notebook index
├── REPOSITORY_GUIDE.md            # This document
├── REVIEW.md                      # Code/science review + fix status
├── query_config.yaml              # Unified config: gene_panel, regions, output, data
├── excitability_genes.md          # Excitability gene mechanisms, tiers, references
├── receptor_excitability.md       # Receptor biology companion (receptors category)
├── pyproject.toml                 # Package metadata (requires-python >=3.11,<3.13)
├── environment.yml                # Conda env (Python 3.12 + pip requirements)
├── requirements.txt               # Core deps
├── .python-version                # pyenv hint (3.12)
│
├── src/                           # Shared library code
│   ├── config.py                  # YAML load/validate, run dirs, parquet discovery
│   ├── data_loaders.py            # ABC Atlas I/O, aggregation, label transfer
│   ├── plotting.py                # Heatmaps, cross-ref scatters, spatial helpers
│   ├── synthesis.py               # Cross-dataset evidence table + confidence tiers
│   └── utils.py                   # Gene ID resolution, brain-area mapping
│
├── notebooks/
│   ├── 01_scrna_heatmaps.ipynb    # M1 — Allen scRNA heatmaps
│   ├── 02_merfish_spatial.ipynb   # M2 — Allen MERFISH heatmaps (filename is legacy)
│   ├── 03_vizgen_crossref.ipynb   # M3 — Vizgen vs Allen cross-ref
│   ├── 04_zhuang_crossref.ipynb   # M4 — Zhuang vs Allen cross-ref
│   └── 05_synthesis.ipynb         # M5 — cross-dataset synthesis
│
├── scripts/
│   └── verify_setup.py            # Smoke test (config, cache, optional Drd2 load)
│
├── data/
│   └── .gitkeep                   # Repo-local cache placeholder only
│
└── figures/                       # Default figure dir (gitignored *.png)
```

**Git ignores:** `figures/*.png`, standard Python/Jupyter artifacts. Downloaded caches live under `/Users/rancze/Documents/Data/expresso_data/`; notebook outputs under `output.output_dir` (see [Outputs and caching](#outputs-and-caching)).

---

## Architecture overview

```
                         query_config.yaml
                    (gene_panel: receptors + excitability)
                              │
                              ▼
                       src/config.py
              (load, validate, derive _all_genes, _categories,
               start_run → timestamped run dir + run_manifest.json)
                              │
     ┌────────────┬───────────┼───────────┬────────────┐
     ▼            ▼           ▼           ▼            ▼
  nb 01        nb 02       nb 03       nb 04        nb 05
  scRNA        MERFISH     Vizgen      Zhuang       synthesis
     │            │           │           │            │
     ▼            ▼           ▼           ▼            ▼
 data_loaders  data_loaders data_loaders data_loaders synthesis.py
     │            │           │           │            │
     ▼            ▼           ▼           ▼            ▼
 plotting.py   plotting.py  plotting.py plotting.py  dot plot +
 (heatmaps)    (heatmaps +   (heatmaps +  (heatmaps +  statement
                cross-ref)    cross-ref)   cross-ref)   scaffold
     │            │           │           │
     └────────────┴───────────┴───────────┘
                         │
              aggregated_*.parquet per run
              (discovered by find_prior_run_parquet)
```

**Design principles:**

| Principle | How it is applied |
|-----------|-------------------|
| Config-driven | Nested `gene_panel` (category → family → genes); brain areas as CCF acronyms; cell type level selectable |
| Partial reads | `anndata.read_h5ad(..., backed='r')` + slice by cells/genes; close file handles in `finally` |
| Graceful degradation | Missing genes warn and skip; empty heatmaps warn and skip |
| Separation of concerns | Notebooks orchestrate; `src/` holds reusable logic |
| Reproducibility | Per-run `run_manifest.json` (git commit + full config snapshot) |
| Cross-dataset honesty | Measured vs imputed provenance tracked; Spearman-first cross-ref; independent-dataset concordance in synthesis |

---

## Configuration reference

File: [`query_config.yaml`](query_config.yaml)

### Required top-level keys

Validated by `src/config.py`:

| Key | Purpose |
|-----|---------|
| `gene_panel` | Nested `category → family → [gene symbols]` (canonical). Legacy top-level `receptors` / `excitability` still accepted. |
| `brain_areas` | List of CCF v3 acronyms |
| `cell_type_level` | One of `class`, `subclass`, `supertype`, `cluster` (default: `supertype`) |
| `output` | Figure settings, output root, caching flags |
| `data` | Cache path, expression unit, MERFISH/Vizgen/Zhuang options |

### Optional top-level keys

| Key | Purpose |
|-----|---------|
| `cell_type_name_filter` | Substrings; keep cell types whose name contains any match (e.g. `["L2/3", "L5"]`) |
| `synthesis` | Optional overrides for detection thresholds (`min_frac`, `min_mean`) in M5 |

### Derived keys (set at load time)

| Key | Type | Description |
|-----|------|-------------|
| `_gene_panel` | `dict[str, list[str]]` | Flattened `family → [genes]` |
| `_gene_panel_key` | `str` | Which top-level key held the panel (`gene_panel`, `receptors`, or `excitability`) |
| `_genes_flat` | `dict[str, str]` | Gene symbol → family name |
| `_all_genes` | `list[str]` | Ordered gene list |
| `_families` | `list[str]` | Family names (unique across categories) |
| `_categories` | `list[str]` | Top-level categories when nested (e.g. `receptors`, `excitability`) |
| `_gene_category` | `dict[str, str]` | Gene symbol → category |
| `_family_category` | `dict[str, str]` | Family → category |
| `_config_path` | `str` | Absolute path to YAML |

`load_config` also accepts a **flat** panel (`family → [genes]`) without categories.

### Current gene panel (committed config)

The unified config defines **243 genes** across **45 families** in **2 categories**:

| Category | Families | Active genes |
|----------|----------|--------------|
| `receptors` | 32 | 146 |
| `excitability` | 13 | 97 |

Many additional genes are listed but commented out in the YAML (e.g. ionotropic glutamate/GABA receptor subunits). Uncomment to activate. To analyse one category only, comment out the other under `gene_panel`.

Per-family heatmaps use the `family` key; synthesis carries the `category` column for receptor vs excitability grouping.

### Brain areas

The config supports two modes (documented inline in the YAML):

1. **Broad divisions** — e.g. `CTX`, `STR`, `TH`, `HY`, `MB`, `HIP`, `AMY`, `CB`, `MY`
2. **CCF sub-regions** — e.g. `VISp`, `VISpm`, `CP`, `CA1`, `DG`

Current setting uses visual cortex sub-regions: `VISp`, `VISpm`, `VISam`, `RSPagl`.

### Output settings

```yaml
output:
  output_dir: /Users/rancze/Documents/!Projects/Ach_NE_Marius_Felix/exploration
  dpi: 150
  heatmap_cmap: viridis
  spatial_cmap: magma
  figsize_heatmap: [14, 8]
  figsize_spatial: [10, 10]
  save_processed_data: true          # Writes aggregated_*.parquet per run
```

Notebooks set `EXPLORATION_ROOT = resolve_output_dir(cfg=config)` after loading the YAML, so `output.output_dir` is the single source of truth. Fallback (if YAML omits it): `DEFAULT_OUTPUT_DIR` in `src/config.py` (currently the `Ach_NE_Marius_Felix/exploration` path).

### Data settings

```yaml
data:
  cache_dir: /Users/rancze/Documents/Data/expresso_data/abc_atlas_cache
  use_imputed_merfish: true          # Fall back to ~8k-gene imputed matrix
  expression_unit: log2                # log2 | raw  → *-log2.h5ad vs *-raw.h5ad
  merfish_dataset: MERFISH-C57BL6J-638850
  vizgen_data_dir: /Users/rancze/Documents/Data/expresso_data/vizgen_cache
  vizgen_samples: null                 # null = all S*R* pairs found in vizgen_data_dir
  vizgen_label_transfer_max_cells: 50000
  vizgen_label_transfer_k: 15
  vizgen_label_transfer_min_confidence: 0.0   # e.g. 0.6 drops low-confidence kNN labels
  zhuang_datasets:                     # Zhuang-ABCA-1 … -4
    - Zhuang-ABCA-1
    - Zhuang-ABCA-2
    - Zhuang-ABCA-3
    - Zhuang-ABCA-4
  allen_merfish_parquet: null          # optional override; auto-discovered if null
```

---

## Data sources

All primary data comes from the ABC Atlas public S3 bucket (`arn:aws:s3:::allen-brain-cell-atlas`), accessed via [`abc_atlas_access`](https://github.com/alleninstitute/abc_atlas_access) (`AbcProjectCache`). Vizgen CSVs are downloaded separately and pointed to by `data.vizgen_data_dir`.

### Datasets used in code

| Directory / source | Dataset | Used for |
|--------------------|---------|----------|
| `WMB-10X` | Metadata hub | Gene table, cell metadata, ROI metadata |
| `WMB-10Xv3` | scRNA-seq expression | M1 (split by anatomical package) |
| `WMB-taxonomy` | Cell type annotations | Join cluster → class/subclass/supertype/cluster |
| `MERFISH-C57BL6J-638850` | MERFISH (~500 genes) | M2 measured expression |
| `MERFISH-C57BL6J-638850-imputed` | Imputed MERFISH (~8k genes) | M2/M3/M4 fallback |
| `MERFISH-C57BL6J-638850-CCF` | CCF coordinates | Spatial helpers in `plotting.py` |
| `Zhuang-ABCA-1` … `-4` | Zhuang MERFISH (~1122 genes) | M4 cross-reference |
| Vizgen CSV pairs | Mouse Brain Receptor Map (~483 genes) | M3 cross-reference |

### Scale (approximate)

| Resource | Size / count | Notes |
|----------|--------------|-------|
| WMB-10Xv3 full brain | ~4M cells, ~32k genes | Split into ~12 anatomical h5ad packages |
| Isocortex packages (example) | ~11.8 GB + ~8.4 GB | Downloaded when querying VIS/CTX cells |
| MERFISH measured matrix | ~7.6 GB | Single h5ad for all cells |
| MERFISH imputed matrix | ~50 GB | Only if genes missing from 500-gene panel |
| MERFISH cells (with CCF) | ~3.74M | After inner join with CCF coords |

Manifest in use (from notebook runs): `releases/20260415/manifest.json`.

### Expression units

Default: **log2(CPM+1)**. Config key `data.expression_unit` selects file suffix (`log2` or `raw`). Plot colorbars label this as `log2(CPM+1)`.

**Vizgen normalisation:** CPM is computed against the **full real-gene panel** (all non-`Blank*` columns), not the requested gene subset. This is critical for cross-dataset comparability (see [REVIEW.md](REVIEW.md) §2.1).

---

## Module reference

### `src/config.py`

| Function | Description |
|----------|-------------|
| `load_config(path)` | Load YAML, validate keys, parse flat or nested `gene_panel`, derive gene lists and category maps |
| `restrict_config_to_genes(config, gene_symbols)` | Prune derived gene lists after a dataset load; returns removed symbols |
| `start_run(project_root, cfg, dataset=..., exploration_root=...)` | Create timestamped run dir + `run_manifest.json` (git + config snapshot) |
| `resolve_output_dir(output_dir, cfg, base_dir)` | Resolve notebook output root (explicit → YAML → `DEFAULT_OUTPUT_DIR`) |
| `get_figures_dir(cfg, base_dir, output_dir)` | Resolve/create figures directory under output root |
| `get_parquet_path(cfg, output_dir, base_dir)` | Path for `aggregated_scrna.parquet` under output root |
| `get_cache_dir(cfg)` | Expands `data.cache_dir` |
| `get_expression_suffix(cfg)` | Returns `'log2'` or `'raw'` |
| `get_vizgen_data_dir(cfg)` | Validates and returns Vizgen CSV directory |
| `get_vizgen_samples(cfg)` | Resolves `vizgen_samples` list (explicit or auto-discovered) |
| `get_zhuang_datasets(cfg)` | Returns configured Zhuang replicate IDs |
| `find_prior_run_parquet(cfg, parquet_filename, dataset_slug, exploration_root)` | Find newest matching aggregate parquet; **warns on `brain_areas` mismatch** vs current config |
| `collect_git_info(project_root)` | Git metadata for run manifests |

Run folder naming: `{timestamp}_{cell_type_level}_{dataset_slug}/`. The manifest stores the config under the legacy key `receptor_query_config` (kept for backward compatibility with existing runs).

### `src/utils.py`

| Function | Description |
|----------|-------------|
| `resolve_gene_ids(gene_df, symbols)` | Map gene symbols → Ensembl IDs; warn on duplicates |
| `warn_missing_genes(found, requested)` | UserWarning for genes not in dataset |
| `build_brain_area_mapping(cache, brain_areas)` | ROI/package → config brain area; returns mapping + assign function |
| `scrna_heatmap_columns(config)` | Column labels for scRNA heatmaps (handles pooled VIS columns) |
| `filter_cell_types_by_name(cell_types, config)` | Apply `cell_type_name_filter` |
| `top_variable_cell_types(matrix, n=50)` | Top-N cell types by row variance (library helper) |
| `assign_merfish_brain_area(...)` | CCF parcellation → config brain area for MERFISH cells |

### `src/data_loaders.py`

**Allen scRNA (M1)**

| Function | Description |
|----------|-------------|
| `get_abc_cache(config)` | `AbcProjectCache.from_cache_dir(...)` |
| `load_scrna_cell_metadata(cache, config)` | WMB-10Xv3 cells + taxonomy + brain_area; filtered to config regions |
| `load_expression_subset(cache, genes, cell_meta, config)` | Partial h5ad load across relevant packages |
| `aggregate_scrna_expression(adata, cell_meta, config)` | Long table: cell_type, brain_area, gene, mean_expression, **frac_expressing**, **n_cells**, family |
| `family_gene_region_matrix(agg_long, family, brain_areas)` | Wide matrix for family heatmap |
| `combined_heatmap_matrix(agg_long, config)` | Wide matrix for combined heatmap |

**Allen MERFISH (M2, M3 reference)**

| Function | Description |
|----------|-------------|
| `load_merfish_cell_metadata(cache, config)` | MERFISH metadata + CCF coords + brain_area |
| `check_gene_availability(cache, gene, config)` | `'present'` \| `'imputed'` \| `'missing'` |
| `check_merfish_genes(cache, genes, config)` | Batch availability check |
| `load_merfish_expression_subset(cache, genes, cell_meta, config)` | Partial multi-gene MERFISH load |
| `aggregate_merfish_expression(adata, cell_meta, config)` | Same schema as scRNA aggregation |
| `family_gene_region_matrix_merfish(agg_long, family, brain_areas)` | MERFISH family matrix |
| `load_allen_merfish_aggregate(config, exploration_root)` | Load cached Allen MERFISH parquet for cross-ref |
| `merfish_gene_source_map(cache, genes, config)` | Gene → `measured` / `imputed` / `missing` |
| `load_single_gene_merfish(cache, gene, config)` | One gene as Series (used by spatial helpers) |

**Zhuang (M4)**

| Function | Description |
|----------|-------------|
| `load_zhuang_cell_metadata(cache, config, dataset_id)` | Zhuang cells + CCF parcellation |
| `check_zhuang_genes(cache, genes, dataset_id)` | Panel availability |
| `load_zhuang_expression_subset(...)` | Partial Zhuang load |
| `aggregate_zhuang_replicates_mean(frames, config)` | Mean across replicates; averages `frac_expressing`, sums `n_cells` |
| `load_zhuang_aggregated(config, exploration_root)` | Full M4 aggregate pipeline |

**Vizgen (M3)**

| Function | Description |
|----------|-------------|
| `load_vizgen_sample(config, sample_tag, genes)` | Load one Vizgen CSV; **full-panel CPM** normalisation |
| `check_vizgen_genes(config, genes)` | Panel availability |
| `build_allen_merfish_label_reference(...)` | Subsampled Allen MERFISH reference for kNN |
| `transfer_allen_merfish_labels(adata, x_ref, y_type, y_area, gene_order, config)` | **Z-scored** kNN label transfer + confidence columns |
| `load_vizgen_aggregated(config, exploration_root)` | Full M3 aggregate pipeline (optional confidence filter) |

**Cross-reference helpers**

| Function | Description |
|----------|-------------|
| `merge_crossref_aggregates(allen_df, other_df, config)` | Inner-join Allen vs other dataset on cell_type × brain_area × gene |

### `src/plotting.py`

| Function | Description |
|----------|-------------|
| `plot_heatmap(agg_df, title, config, save_path, ...)` | Seaborn clustermap (or plain heatmap if <2×2) |
| `plot_family_heatmap(family, gene_matrix, config, ...)` | Per-family heatmap |
| `plot_combined_heatmap(all_genes_matrix, config, ...)` | Combined all-genes heatmap |
| `plot_crossref_scatter(merged, config, title, save_path, ...)` | **Spearman ρ-first** Allen vs other scatter |
| `plot_crossref_family_scatters(...)` | Per-family cross-ref scatters |
| `plot_crossref_side_by_side_heatmaps(...)` | Aligned Allen vs other heatmaps |
| `plot_spatial(coords_df, expression, gene, projection, config, ...)` | Single-gene CCF scatter (library extra; not wired into NB02) |
| `plot_family_spatial_panel(...)` | Grid: genes × projections (library extra) |
| `format_gene_label(gene, source)` | Appends `*` for imputed genes |

### `src/synthesis.py`

| Function | Description |
|----------|-------------|
| `discover_dataset_parquets(config, exploration_root)` | Locate newest `aggregated_*.parquet` per dataset |
| `gather_dataset_aggregates(config, exploration_root)` | Read all available aggregates; warn on missing |
| `build_evidence_table(aggregates, config, allen_gene_sources=...)` | Core evidence table: (cell_type, brain_area, gene) with per-dataset metrics, detection flags, confidence tier |
| `target_evidence(evidence, cell_type, region=None)` | Filter to one cell type (± region) |
| `summarize_target(evidence, cell_type, region=None)` | Tier/family/category summary for a target |
| `statement_scaffold(summary)` | Human-readable claim split by category |
| `plot_evidence_dotplot(evidence, target, config, save_path)` | Dot plot: mean × fraction across datasets |

`DATASET_SPECS` registry defines per-dataset parquet names, slug matching, region resolution, and imputed support. Allen MERFISH imputed values are **supporting** evidence only (not independent of Allen scRNA).

---

## Milestone 1 — scRNA-seq heatmaps

**Notebook:** [`notebooks/01_scrna_heatmaps.ipynb`](notebooks/01_scrna_heatmaps.ipynb)

### Pipeline steps

1. **Load config** → derive gene list, brain areas, cell type level
2. **`start_run`** → timestamped output dir + manifest
3. **Init cache** → `get_abc_cache(config)`
4. **Load cell metadata** → `load_scrna_cell_metadata`
   - Filter to `dataset_label == "WMB-10Xv3"`
   - Join WMB taxonomy on `cluster_alias`
   - Assign `brain_area` via ROI/package mapping
   - Filter to config `brain_areas`
5. **Load expression** → `load_expression_subset` (backed partial reads per package)
6. **Aggregate** → `aggregate_scrna_expression`
   - Mean expression, fraction expressing, and cell count per `(cell_type_level, brain_area, gene)`
7. **Plot** per-family and combined heatmaps
8. **Cache** → `aggregated_scrna.parquet` in the run dir

### Aggregated output schema

| Column | Description |
|--------|-------------|
| `cell_type` | Value at configured taxonomy level (e.g. supertype name) |
| `brain_area` | Config brain area acronym (may be pooled — see §Brain area mapping) |
| `gene` | Gene symbol |
| `mean_expression` | Mean log2(CPM+1) across cells in group (zeros included) |
| `frac_expressing` | Fraction of cells with expression > 0 |
| `n_cells` | Number of cells in the group |
| `family` | Gene family from config |

### Partial loading pattern (critical)

The full WMB-10Xv3 matrix is far too large for RAM. The code filters cells first (metadata only), loads only packages present in `cell_meta.feature_matrix_label`, uses backed mode, and closes files after each package.

---

## Brain area mapping (scRNA caveat)

scRNA-seq cells in WMB-10X are annotated with **dissection ROI** (`region_of_interest_acronym`), not per-cell CCF parcellation. Expresso bridges this gap in `build_brain_area_mapping`:

```
Config brain_area (CCF-style)     scRNA dissection ROI / package
─────────────────────────────     ───────────────────────────────
VISp, VISpm, VISam, ...      →    VIS  (with UserWarning)
CP, ACB                      →    STRd, STRv → STR
CA1, DG                      →    HIP
CTX (broad)                  →    All cortical ROIs in _CORTICAL_ROIS
TH, HY, MB, ...              →    Direct acronym match or package suffix
```

**Implications:**

- Querying `VISp` alone still pulls **all VIS dissection cells**, then labels them as `VISp`. You cannot resolve primary vs secondary visual cortex at scRNA resolution.
- Warnings are intentional — read them when interpreting heatmaps.
- MERFISH (M2), Vizgen (M3), and Zhuang (M4) **do** have CCF parcellation and support fine-region claims.
- Synthesis flags Allen scRNA as `region_resolved=False` and broadcasts its values across regions.

---

## Milestone 2 — MERFISH heatmaps

**Notebook:** [`notebooks/02_merfish_spatial.ipynb`](notebooks/02_merfish_spatial.ipynb) (filename reflects an earlier spatial-scatter design; the notebook now produces **heatmaps** like M1 but with native CCF sub-region resolution).

### Pipeline steps

1. Load config, `start_run`, init cache
2. **Load MERFISH metadata** → `load_merfish_cell_metadata` (CCF parcellation → `brain_area`)
3. **Check genes** → `check_merfish_genes` (panel / imputed / missing)
4. **Load expression** → `load_merfish_expression_subset` (batched partial read)
5. **Aggregate** → `aggregate_merfish_expression` (same schema as M1, including `frac_expressing` / `n_cells`)
6. **Plot** per-family and combined heatmaps
7. **Cache** → `aggregated_merfish.parquet`

### Gene availability logic

```
Gene in 500-gene panel?  →  load from MERFISH-C57BL6J-638850  (measured)
Else if use_imputed_merfish and in imputed var?  →  load from -imputed  (imputed)
Else  →  skip with warning
```

Imputed genes are marked with `*` in heatmap labels. Most excitability ion-channel genes require the imputed matrix or are absent from MERFISH entirely.

`plot_spatial` / `plot_family_spatial_panel` remain in `plotting.py` as optional library helpers (single-gene CCF scatter maps) but are not called by the current notebook.

---

## Milestones 3–4 — Cross-reference

### Milestone 3 — Vizgen (`03_vizgen_crossref.ipynb`)

- **Data:** Vizgen MERFISH Mouse Brain Receptor Map (~483 genes); flat CSV pairs under `data.vizgen_data_dir`
- **Normalisation:** full real-gene panel CPM (exclude `Blank*` probes)
- **Label transfer:** z-scored kNN on shared genes transfers Allen `supertype` and CCF `brain_area` to Vizgen cells; confidence recorded per cell; optional filter via `vizgen_label_transfer_min_confidence`
- **Outputs:** Vizgen heatmaps, Allen↔Vizgen scatter plots (Spearman ρ-first), side-by-side heatmaps, `aggregated_vizgen.parquet`
- **Coverage note:** Vizgen is a receptor panel — almost no excitability genes

### Milestone 4 — Zhuang (`04_zhuang_crossref.ipynb`)

- **Data:** `Zhuang-ABCA-1` … `Zhuang-ABCA-4` via `AbcProjectCache`; ~1122-gene panel
- **Aggregation:** replicate mean with averaged `frac_expressing` and summed `n_cells`
- **Outputs:** Zhuang heatmaps, Allen↔Zhuang scatter plots, `aggregated_zhuang.parquet`
- **Reuse:** discovers cached `aggregated_merfish.parquet` via `find_prior_run_parquet` (warns on `brain_areas` mismatch)

Cross-ref plots lead with **Spearman ρ** (rank concordance); Pearson r is secondary. Measured-only ρ is reported separately from imputed Allen genes.

---

## Milestone 5 — Synthesis

**Notebook:** [`notebooks/05_synthesis.ipynb`](notebooks/05_synthesis.ipynb)  
**Module:** [`src/synthesis.py`](src/synthesis.py)

### Pipeline steps

1. Load unified `query_config.yaml` (receptors + excitability)
2. **`gather_dataset_aggregates`** — discover and read `aggregated_scrna/merfish/vizgen/zhuang.parquet` from prior runs (warns if `brain_areas` differ)
3. **`merfish_gene_source_map`** — Allen MERFISH measured vs imputed provenance
4. **`build_evidence_table`** — one row per (cell_type, brain_area, gene) with:
   - per-dataset `mean_expression`, `frac_expressing`, `n_cells`, `source` (measured/imputed)
   - detection flag (`frac ≥ min_frac` AND `mean ≥ min_mean`; defaults 0.25 / 1.0, overridable via optional `synthesis:` config block)
   - `n_independent_measured_detections` and `confidence_tier` (high / medium / low / none)
   - `category` column (receptors vs excitability)
5. **Target analysis** — `target_evidence` / `summarize_target` / `statement_scaffold` for a chosen cell type (± region)
6. **Plot** → `plot_evidence_dotplot`
7. **Write** → `evidence_table.parquet` + `.csv`

Level-agnostic: works at whatever `cell_type_level` is configured. Re-run notebooks 01–04 after code changes so parquets carry `frac_expressing` and `n_cells`.

---

## Plotting layer

### Heatmaps

- Uses `seaborn.clustermap` when matrix is ≥2×2 (row and column clustering)
- Falls back to plain `sns.heatmap` for tiny matrices
- Colormap from `output.heatmap_cmap` (default `viridis`)
- Titles are panel-neutral ("Family — mean log2(CPM+1)", not hardcoded "receptors")
- Imputed genes tagged with `*` via `format_gene_label`

### Cross-reference

- `plot_crossref_scatter`: Spearman ρ in title (overall + measured-only); Pearson r secondary
- Hollow markers / `*` for imputed Allen genes
- Side-by-side heatmaps align cell types and genes between datasets

### Spatial (library extras)

- `plot_spatial` / `plot_family_spatial_panel`: CCF scatter maps; not wired into current notebooks
- 99th percentile of positive values for `vmax`; imputed genes tagged in title

---

## Notebooks

| Notebook | Status | Primary outputs |
|----------|--------|-----------------|
| `01_scrna_heatmaps.ipynb` | ✅ Implemented | `heatmap_{family}.png`, `heatmap_combined.png`, `aggregated_scrna.parquet` |
| `02_merfish_spatial.ipynb` | ✅ Implemented | `heatmap_{family}.png`, `heatmap_combined.png`, `aggregated_merfish.parquet` |
| `03_vizgen_crossref.ipynb` | ✅ Implemented | Vizgen heatmaps + Allen cross-ref scatters/heatmaps, `aggregated_vizgen.parquet` |
| `04_zhuang_crossref.ipynb` | ✅ Implemented | Zhuang heatmaps + Allen cross-ref scatters, `aggregated_zhuang.parquet` |
| `05_synthesis.ipynb` | ✅ Implemented | `evidence_table.parquet`/`.csv`, focused dot plot, statement scaffold |

### Common notebook boilerplate

All notebooks resolve `PROJECT_ROOT` (cwd or parent if run from `notebooks/`), insert into `sys.path`, load `query_config.yaml`, set `EXPLORATION_ROOT = resolve_output_dir(cfg=config)`, call `start_run(..., exploration_root=EXPLORATION_ROOT)`, and import from `src.*`. Kernelspec: conda env `expresso`.

### Dev / smoke cell (M1 only)

The last cells in `01_scrna_heatmaps.ipynb` override config to **Drd1/Drd2 × STR/TH** for a smaller download footprint.

---

## Outputs and caching

| Location | Config key | Contents |
|----------|------------|----------|
| `…/Ach_NE_Marius_Felix/exploration` | `output.output_dir` | Timestamped run dirs, parquets, figures |
| `…/Data/expresso_data/abc_atlas_cache` | `data.cache_dir` | ABC Atlas downloads (scRNA, MERFISH, Zhuang) |
| `…/Data/expresso_data/vizgen_cache` | `data.vizgen_data_dir` | Vizgen MERFISH CSV pairs |

| Artifact | Location | When |
|----------|----------|------|
| ABC Atlas files | `data.cache_dir` | First access per file |
| Run directory | `{exploration_root}/{timestamp}_{level}_{dataset}/` | Each notebook run via `start_run` |
| Run manifest | `{run_dir}/run_manifest.json` | Git commit + full config snapshot |
| Aggregated tables | `{run_dir}/aggregated_{scrna,merfish,vizgen,zhuang}.parquet` | M1–M4 |
| Evidence table | `{run_dir}/evidence_table.parquet` + `.csv` | M5 |
| Heatmaps / scatters | `{run_dir}/` or figures subdir | M1–M5 |
| Verify test heatmap | `figures/verify_test_heatmap.png` | `verify_setup.py` |

`find_prior_run_parquet` discovers the **newest** matching aggregate under the exploration root (matches `cell_type_level` + `dataset_slug`). It issues a **`REGION MISMATCH` warning** when the discovered run's `brain_areas` differ from the current config — cross-ref and synthesis only join on overlapping regions.

Parquet enables re-plotting and synthesis without re-downloading expression matrices. Saved notebook outputs in the repo may reflect an older config or code version; re-execute after YAML or code edits.

---

## Environment and dependencies

| Requirement | Detail |
|-------------|--------|
| Python | **3.11 or 3.12** only (`requires-python = ">=3.11,<3.13"`) |
| Recommended setup | `conda env create -f environment.yml && conda activate expresso` |
| Pip alternative | `pip install -r requirements.txt` |
| Vizgen (M3) | Same as M1/M2; Vizgen CSV I/O uses pandas/anndata |

### Core packages

- `abc_atlas_access` (from GitHub)
- `anndata`, `pandas`, `numpy`, `scipy`
- `matplotlib`, `seaborn`
- `pyyaml`, `pyarrow`, `tqdm`

---

## Smoke testing

```bash
python scripts/verify_setup.py --quick   # Config, cache init, mock heatmap
python scripts/verify_setup.py             # Also loads Drd2 from WMB-10Xv3-STR
```

Steps:

1. Load `query_config.yaml`; assert genes present
2. Init `AbcProjectCache`
3. Check Drd2 in gene metadata
4. (Full) Load STR cells + Drd2 expression via backed read
5. Render 2×2 test heatmap to `figures/verify_test_heatmap.png`

Exits with error on unsupported Python (≥3.13 or <3.11).

---

## Known limitations and gotchas

1. **scRNA sub-regions are approximate** — CCF acronyms in config map to coarser dissection ROIs; warnings indicate when this happens.
2. **Large first-run downloads** — Plan disk space (~20 GB+ for cortex scRNA; ~50 GB if imputed MERFISH needed).
3. **MERFISH panel coverage** — Only ~500 genes measured; many receptors and most excitability genes require imputed matrix or scRNA-only analysis.
4. **Circular validation** — Allen MERFISH imputed values are predicted from Allen scRNA; treat them as supporting, not independent, evidence.
5. **Cross-dataset magnitudes differ** — Panels differ in size; use Spearman ρ and detection concordance, not raw log2(CPM+1) equality across datasets.
6. **Vizgen label transfer is approximate** — Vizgen cells get Allen labels via kNN, not native Vizgen annotations; use `vizgen_label_transfer_min_confidence` to filter.
7. **Region mismatch across runs** — `find_prior_run_parquet` warns but does not block; re-run source notebooks when changing `brain_areas`.
8. **Saved notebook outputs may be stale** — Re-execute after YAML or code edits; re-run M1–M4 before M5 if parquets lack `frac_expressing`/`n_cells`.

---

## Suggested workflows

### Quick exploration (minimal download)

1. Narrow `query_config.yaml` to 1–2 genes and 1–2 broad regions (`STR`, `TH`)
2. Run M1 dev cell or `verify_setup.py` (without `--quick`)
3. Run M2 for genes in the 500-gene panel

### Full receptor + excitability survey

1. Use the full `gene_panel` (or comment out one category to focus)
2. Set `brain_areas` and `cell_type_level` (default `supertype`; optional `cell_type_name_filter`)
3. Run `01_scrna_heatmaps.ipynb` — expect long runtime and large downloads for cortical regions
4. Run `02_merfish_spatial.ipynb` — set `use_imputed_merfish: false` first to avoid 50 GB download; enable for missing genes

### Cross-dataset validation and synthesis

1. Run M1 + M2 (Allen scRNA + MERFISH)
2. Run M3 (Vizgen) and M4 (Zhuang) — they discover cached Allen MERFISH aggregates
3. Run M5 (synthesis) to build the evidence table, concordance, confidence tiers, and statement scaffold for a target cell type (e.g. an L5 supertype in `VISpm`)

---

## Document map

| Document | Audience | Content |
|----------|----------|---------|
| [README.md](README.md) | New users | Setup, notebook list, smoke test |
| **REPOSITORY_GUIDE.md** | Developers / analysts | Architecture, APIs, data flows, caveats |
| [REVIEW.md](REVIEW.md) | Maintainers | Code/science review, fix status |
| [excitability_genes.md](excitability_genes.md) | Analysts | Excitability gene mechanisms, tiers, references |
| [receptor_excitability.md](receptor_excitability.md) | Analysts | Receptor biology companion for the receptors category |

---

*ABC Atlas manifest and download sizes reflect runs against `releases/20260415`.*
