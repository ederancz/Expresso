# Expresso — Repository Deep Dive

A technical guide to the **Expresso** codebase: purpose, architecture, data flows, module APIs, configuration, and operational notes. For quick setup, see [README.md](README.md). For milestone implementation sketches aimed at development agents, see [cursor_handover.md](cursor_handover.md).

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
9. [Milestone 2 — MERFISH spatial maps](#milestone-2--merfish-spatial-maps)
10. [Plotting layer](#plotting-layer)
11. [Notebooks](#notebooks)
12. [Outputs and caching](#outputs-and-caching)
13. [Environment and dependencies](#environment-and-dependencies)
14. [Smoke testing](#smoke-testing)
15. [Future milestones](#future-milestones)
16. [Known limitations and gotchas](#known-limitations-and-gotchas)
17. [Suggested workflows](#suggested-workflows)

---

## What this project does

Expresso queries **receptor gene expression** across **mouse brain areas** and **cell types** using the [Allen Brain Cell Atlas (ABC Atlas)](https://alleninstitute.github.io/abc_atlas_access/). It is designed as a **config-driven analysis pipeline**:

- **No hard-coded genes or regions** in notebook logic — everything comes from `receptor_query_config.yaml`.
- **Memory-efficient I/O** — expression matrices are read in backed mode and sliced to only the genes and cells needed.
- **Two primary analysis modes** (implemented):
  - **Milestone 1:** scRNA-seq heatmaps (cell type × brain area)
  - **Milestone 2:** MERFISH spatial scatter maps (whole-brain CCF coordinates)
- **Two cross-reference modes** (planned):
  - **Milestone 3:** Vizgen MERFISH Mouse Receptor Map
  - **Milestone 4:** Zhuang MERFISH replicates via ABC Atlas

The project name reflects its focus: **express**ion of neurotransmitter **receptors** (dopamine, serotonin, glutamate, GABA, opioid, cannabinoid, acetylcholine, adrenergic families).

---

## Repository layout

```
Expresso/
├── README.md                      # Quick start, setup, notebook index
├── REPOSITORY_GUIDE.md            # This document
├── cursor_handover.md             # Implementation notes for AI / future milestones
├── receptor_query_config.yaml     # Runtime config (genes, regions, output, data paths)
├── pyproject.toml                 # Package metadata (requires-python >=3.11,<3.13)
├── environment.yml                # Conda env (Python 3.12 + pip requirements)
├── requirements.txt               # Core deps (Milestones 1–3)
├── .python-version                # pyenv hint (3.12)
│
├── src/                           # Shared library code
│   ├── config.py                  # YAML load/validate, path helpers
│   ├── data_loaders.py            # ABC Atlas I/O, aggregation
│   ├── plotting.py                # Heatmaps and spatial plots
│   └── utils.py                   # Gene ID resolution, brain-area mapping
│
├── notebooks/
│   ├── 01_scrna_heatmaps.ipynb    # Milestone 1 — implemented, run successfully
│   ├── 02_merfish_spatial.ipynb   # Milestone 2 — implemented, run successfully
│   ├── 03_vizgen_crossref.ipynb   # Milestone 3 — Vizgen cross-ref vs Allen
│   └── 04_zhuang_crossref.ipynb   # Milestone 4 — TODO stub
│
├── scripts/
│   └── verify_setup.py            # Smoke test (config, cache, optional Drd2 load)
│
├── data/
│   ├── .gitkeep
│   └── aggregated_scrna.parquet   # Cached M1 output (when save_processed_data: true)
│
└── figures/                       # Default figure dir (gitignored *.png); notebooks may override
```

**Git ignores:** `figures/*.png`, standard Python/Jupyter artifacts. ABC Atlas cache (`~/abc_atlas_cache` by default) lives outside the repo.

---

## Architecture overview

```
                    receptor_query_config.yaml
                              │
                              ▼
                       src/config.py
                    (load, validate, derive
                     _all_genes, _families)
                              │
              ┌───────────────┴───────────────┐
              ▼                               ▼
     notebooks/01_scrna              notebooks/02_merfish
              │                               │
              ▼                               ▼
     src/data_loaders.py              src/data_loaders.py
     (WMB-10Xv3 partial load)        (single-gene MERFISH load)
              │                               │
              ▼                               ▼
     src/plotting.py                  src/plotting.py
     (heatmaps)                       (spatial scatter + panels)
              │                               │
              ▼                               ▼
     figures/, data/*.parquet          figures/spatial_*.png
```

**Design principles:**

| Principle | How it is applied |
|-----------|-------------------|
| Config-driven | Genes grouped by receptor family; brain areas as CCF acronyms; cell type level selectable |
| Partial reads | `anndata.read_h5ad(..., backed='r')` + slice by cells/genes; close file handles in `finally` |
| Graceful degradation | Missing genes warn and skip; empty heatmaps warn and skip |
| Separation of concerns | Notebooks orchestrate; `src/` holds reusable logic |
| External cache | `AbcProjectCache` downloads to `data.cache_dir` (default `~/abc_atlas_cache`) |

---

## Configuration reference

File: [`receptor_query_config.yaml`](receptor_query_config.yaml)

### Required top-level keys

Validated by `src/config.py`:

| Key | Purpose |
|-----|---------|
| `receptors` | Dict of family name → list of gene symbols |
| `brain_areas` | List of CCF v3 acronyms |
| `cell_type_level` | One of `class`, `subclass`, `supertype`, `cluster` |
| `output` | Figure settings, caching flags |
| `data` | Cache path, expression unit, MERFISH options |

### Derived keys (set at load time)

| Key | Type | Description |
|-----|------|-------------|
| `_genes_flat` | `dict[str, str]` | Gene symbol → family name |
| `_all_genes` | `list[str]` | Ordered gene list |
| `_families` | `list[str]` | Receptor family names |
| `_config_path` | `str` | Absolute path to YAML (for relative path resolution) |

### Current gene panel (full config)

The committed config defines **~50 genes** across **8 families**:

| Family | Genes |
|--------|-------|
| dopamine | Drd1, Drd2, Drd3, Drd4, Drd5 |
| serotonin | Htr1a, Htr1b, Htr2a, Htr2c, Htr3a, Htr4, Htr6, Htr7 |
| glutamate | Grin1, Grin2a, Grin2b, Grm1, Grm5 |
| gaba | Gabra1, Gabra2, Gabrb2, Gabrg2 |
| opioid | Oprm1, Oprd1, Oprk1 |
| cannabinoid | Cnr1, Cnr2 |
| acetylcholine | Chrm1–5, Chrna2–4, Chrna7, Chrnb2, Chrnb4 |
| adrenergic | Adra1a/b/d, Adra2a/b/c, Adrb1–3 |

Note: notebook execution outputs in the repo reflect an **earlier narrowed config** (e.g. only `Htr2a` in `VISp`). Re-run notebooks after editing the YAML to pick up changes.

### Brain areas

The config supports two modes (documented inline in the YAML):

1. **Broad divisions** — e.g. `CTX`, `STR`, `TH`, `HY`, `MB`, `HIP`, `AMY`, `CB`, `MY`
2. **CCF sub-regions** — e.g. `VISp`, `VISpm`, `CP`, `CA1`, `DG`

Current setting uses visual cortex sub-regions: `VISp`, `VISpm`, `VISam`, `RSPagl`.

### Output settings

```yaml
output:
  figures_dir: figures/           # Fallback; notebooks override via OUTPUT_DIR
  dpi: 150
  heatmap_cmap: viridis
  spatial_cmap: magma
  figsize_heatmap: [14, 8]
  figsize_spatial: [10, 10]
  save_processed_data: true        # Writes data/aggregated_scrna.parquet
```

### Data settings

```yaml
data:
  cache_dir: ~/abc_atlas_cache
  use_imputed_merfish: true        # Fall back to ~8k-gene imputed matrix
  expression_unit: log2            # log2 | raw  → selects *-log2.h5ad vs *-raw.h5ad
  merfish_dataset: MERFISH-C57BL6J-638850
  vizgen_data_dir: null            # For Milestone 3
```

---

## Data sources

All primary data comes from the ABC Atlas public S3 bucket (`arn:aws:s3:::allen-brain-cell-atlas`), accessed via [`abc_atlas_access`](https://github.com/alleninstitute/abc_atlas_access) (`AbcProjectCache`).

### Datasets used in code

| Directory constant | Dataset | Used for |
|--------------------|---------|----------|
| `WMB-10X` | Metadata hub | Gene table, cell metadata, ROI metadata |
| `WMB-10Xv3` | scRNA-seq expression | Milestone 1 (split by anatomical package) |
| `WMB-taxonomy` | Cell type annotations | Join cluster → class/subclass/supertype/cluster |
| `MERFISH-C57BL6J-638850` | MERFISH (~500 genes) | Milestone 2 measured expression |
| `MERFISH-C57BL6J-638850-imputed` | Imputed MERFISH (~8k genes) | Milestone 2 fallback |
| `MERFISH-C57BL6J-638850-CCF` | CCF coordinates | Spatial plotting (`x_ccf`, `y_ccf`, `z_ccf`) |

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

---

## Module reference

### `src/config.py`

| Function | Description |
|----------|-------------|
| `load_config(path)` | Load YAML, validate required keys and `cell_type_level`, derive gene lists |
| `get_figures_dir(cfg, base_dir, output_dir)` | Resolve/create figures directory; `output_dir` overrides config |
| `get_parquet_path(cfg, base_dir)` | Returns `{base_dir}/data/aggregated_scrna.parquet` |
| `get_cache_dir(cfg)` | Expands `data.cache_dir` |
| `get_expression_suffix(cfg)` | Returns `'log2'` or `'raw'` |

### `src/utils.py`

| Function | Description |
|----------|-------------|
| `resolve_gene_ids(gene_df, symbols)` | Map gene symbols → Ensembl IDs; warn on duplicates |
| `warn_missing_genes(found, requested)` | UserWarning for genes not in dataset |
| `build_brain_area_mapping(cache, brain_areas)` | ROI/package → config brain area; returns mapping + assign function |
| `top_variable_cell_types(matrix, n=50)` | Top-N cell types by row variance (for combined heatmap) |

**Internal mappings:**

- `_CORTICAL_ROIS` — dissection acronyms grouped as `CTX`
- `_CCF_TO_DISSECTION_ROI` — CCF sub-regions (e.g. `VISp` → `VIS`, `CP` → `STRd`)
- `_PACKAGE_TO_BRAIN_AREA` — `feature_matrix_label` suffix → brain area (e.g. `Isocortex-1` → `CTX`)

### `src/data_loaders.py`

| Function | Description |
|----------|-------------|
| `get_abc_cache(config)` | `AbcProjectCache.from_cache_dir(...)` |
| `load_scrna_cell_metadata(cache, config)` | WMB-10Xv3 cells + taxonomy + brain_area; filtered to config regions |
| `load_expression_subset(cache, genes, cell_meta, config)` | Partial h5ad load across relevant packages |
| `aggregate_scrna_expression(adata, cell_meta, config)` | Long DataFrame: cell_type, brain_area, gene, mean_expression, family |
| `family_gene_region_matrix(agg_long, family, brain_areas)` | Wide matrix for family heatmap (cell types × brain areas) |
| `combined_heatmap_matrix(agg_long, top_cell_types, brain_areas)` | Wide matrix for combined heatmap (cell types × genes) |
| `load_merfish_cell_metadata(cache, config)` | MERFISH metadata + CCF coords |
| `check_gene_availability(cache, gene, config)` | `'present'` \| `'imputed'` \| `'missing'` |
| `load_single_gene_merfish(cache, gene, config)` | One gene as Series; source `'measured'` or `'imputed'` |

### `src/plotting.py`

| Function | Description |
|----------|-------------|
| `plot_heatmap(agg_df, title, config, save_path, ...)` | Seaborn clustermap (or plain heatmap if <2×2) |
| `plot_family_heatmap(family, gene_matrix, config, ...)` | Saves `heatmap_{family}.png` |
| `plot_combined_heatmap(all_genes_matrix, config, ...)` | Saves `heatmap_combined.png` |
| `plot_spatial(coords_df, expression, gene, projection, config, ...)` | Single-gene CCF or section scatter |
| `plot_family_spatial_panel(family_results, family, coords_df, config, ...)` | Grid: genes × projections |

---

## Milestone 1 — scRNA-seq heatmaps

**Notebook:** [`notebooks/01_scrna_heatmaps.ipynb`](notebooks/01_scrna_heatmaps.ipynb)

### Pipeline steps

1. **Load config** → derive gene list, brain areas, cell type level
2. **Init cache** → `get_abc_cache(config)`
3. **Load cell metadata** → `load_scrna_cell_metadata`
   - Filter to `dataset_label == "WMB-10Xv3"`
   - Join WMB taxonomy on `cluster_alias`
   - Assign `brain_area` via ROI/package mapping
   - Filter to config `brain_areas`
4. **Load expression** → `load_expression_subset`
   - Resolve symbols → Ensembl IDs via WMB-10X gene table
   - Iterate unique `feature_matrix_label` values in filtered cells
   - For each package: backed read `{pkg}/{log2|raw}`, slice cells × genes, `to_memory()`
   - Concatenate packages; attach `gene_symbol` to `var`
5. **Aggregate** → `aggregate_scrna_expression`
   - Mean expression per `(cell_type_level, brain_area)` for each gene
6. **Plot**
   - Per family: `family_gene_region_matrix` → `plot_family_heatmap`
   - Combined: top 50 variable cell types → `plot_combined_heatmap`
7. **Cache** (optional) → `agg_long.to_parquet(data/aggregated_scrna.parquet)`

### Aggregated output schema

| Column | Description |
|--------|-------------|
| `cell_type` | Value at configured taxonomy level (e.g. supertype name) |
| `brain_area` | Config brain area acronym |
| `gene` | Gene symbol |
| `mean_expression` | Mean log2(CPM+1) across cells in group |
| `family` | Receptor family from config |

Example run (Htr2a, VISp): **30,882 cells**, **102 aggregated rows** (102 supertypes × 1 gene × 1 region).

### Partial loading pattern (critical)

The full WMB-10Xv3 matrix is far too large for RAM. The code:

- Filters cells **first** (metadata only)
- Loads only **packages** present in `cell_meta.feature_matrix_label`
- Uses **backed mode** and closes files after each package
- Converts only the sliced subset to memory

This mirrors the pattern documented in `cursor_handover.md` and Allen Institute tutorials.

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

- Querying `VISp` alone still pulls **all VIS dissection cells** (~31k), then labels them as `VISp`. You cannot resolve primary vs secondary visual cortex at scRNA resolution without MERFISH/spatial data.
- Warnings are intentional — read them when interpreting heatmaps.
- MERFISH (Milestone 2) **does** have CCF coordinates and is the right modality for sub-regional spatial patterns.

---

## Milestone 2 — MERFISH spatial maps

**Notebook:** [`notebooks/02_merfish_spatial.ipynb`](notebooks/02_merfish_spatial.ipynb)

### Pipeline steps

1. Load config and cache
2. **Load MERFISH metadata** → `load_merfish_cell_metadata`
   - `cell_metadata_with_cluster_annotation` from MERFISH dataset
   - Join `ccf_coordinates` from `MERFISH-C57BL6J-638850-CCF`
   - Rename section coords to `x_section`, `y_section`, `z_section`
   - CCF coords: `x_ccf`, `y_ccf`, `z_ccf`
3. For each gene in config:
   - `check_gene_availability` → panel / imputed / missing
   - `load_single_gene_merfish` → expression Series (~3.74M values)
   - Plot coronal, sagittal, axial via `plot_spatial(coord_prefix='ccf')`
4. **Family panels** → `plot_family_spatial_panel` (genes × 3 projections)

### Gene availability logic

```
Gene in 500-gene panel?  →  load from MERFISH-C57BL6J-638850  (measured)
Else if use_imputed_merfish and in imputed var?  →  load from -imputed  (imputed)
Else  →  skip with warning
```

Imputed gene symbols are cached module-wide in `_imputed_panel_cache` after first backed read of imputed h5ad `var`.

**Note:** `load_single_gene_merfish` reads the **full MERFISH h5ad** per gene (backed, one column). For many genes this is I/O-heavy but RAM-light. A future optimization could batch genes or memory-map more aggressively.

### Projections

| Projection | CCF axes | Section axes (unused by default) |
|------------|----------|----------------------------------|
| coronal | x_ccf, y_ccf | x_section, y_section |
| sagittal | z_ccf, y_ccf | z_section, y_section |
| axial | x_ccf, z_ccf | x_section, z_section |

Plots use 99th percentile of positive values for `vmax`, `vmin=0`, rasterized scatter, equal aspect, inverted y-axis.

---

## Plotting layer

### Heatmaps

- Uses `seaborn.clustermap` when matrix is ≥2×2 (row and column clustering)
- Falls back to plain `sns.heatmap` for tiny matrices
- Colormap from `output.heatmap_cmap` (default `viridis`)

### Spatial

- Point size 0.3 (single) / 0.2 (panel), alpha 0.6 / 0.5
- Imputed genes tagged in title: `[imputed]` or `[imp]` in panels
- Family panel figure size scales with gene count and projection count

### Figure output paths

Notebooks set:

```python
OUTPUT_DIR = Path("/Users/rancze/Documents/!Projects/Ach_NE_Marius_Felix/exploration")
figures_dir = get_figures_dir(config, output_dir=OUTPUT_DIR)
```

Some notebook cells print paths under `Expresso/figures/` — that happens when `get_figures_dir` resolves differently or outputs were from an earlier run. **Edit `OUTPUT_DIR` in each notebook** to control where PNGs land.

---

## Notebooks

| Notebook | Status | Primary outputs |
|----------|--------|-----------------|
| `01_scrna_heatmaps.ipynb` | ✅ Complete | `heatmap_{family}.png`, `heatmap_combined.png`, optional parquet |
| `02_merfish_spatial.ipynb` | ✅ Complete | `spatial_{gene}_{projection}.png`, `spatial_panel_{family}.png` |
| `03_vizgen_crossref.ipynb` | ✅ Implemented | Vizgen heatmaps + Allen cross-ref (parquet reuse) |
| `04_zhuang_crossref.ipynb` | 🚧 Stub | Config load + OUTPUT_DIR only |

### Common notebook boilerplate

All notebooks resolve `PROJECT_ROOT` (cwd or parent if run from `notebooks/`), insert into `sys.path`, and import from `src.*`. Kernelspec: conda env `expresso`.

### Dev / smoke cell (M1 only)

The last cells in `01_scrna_heatmaps.ipynb` override config to **Drd1/Drd2 × STR/TH** for a smaller download footprint — useful for validation without pulling full isocortex packages.

---

## Outputs and caching

| Artifact | Location | When |
|----------|----------|------|
| ABC Atlas files | `~/abc_atlas_cache` (configurable) | First access per file |
| Aggregated scRNA table | `data/aggregated_scrna.parquet` | M1, if `save_processed_data: true` |
| Heatmaps | `OUTPUT_DIR` or `figures/` | M1 |
| Spatial PNGs | `OUTPUT_DIR` or `figures/` | M2 |
| Verify test heatmap | `figures/verify_test_heatmap.png` | `verify_setup.py` |

Parquet enables re-plotting or downstream analysis without re-downloading expression matrices.

---

## Environment and dependencies

| Requirement | Detail |
|-------------|--------|
| Python | **3.11 or 3.12** only (`requires-python = ">=3.11,<3.13"`) |
| Recommended setup | `conda env create -f environment.yml && conda activate expresso` |
| Pip alternative | `pip install -r requirements.txt` |
| Vizgen (M3) | Same as M1/M2 (`requirements.txt`); Vizgen CSV I/O uses pandas/anndata |

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

1. Load config; assert genes present
2. Init `AbcProjectCache`
3. Check Drd2 in gene metadata
4. (Full) Load STR cells + Drd2 expression via backed read
5. Render 2×2 test heatmap to `figures/verify_test_heatmap.png`

Exits with error on unsupported Python (≥3.13 or <3.11).

---

## Future milestones

### Milestone 3 — Vizgen (`03_vizgen_crossref.ipynb`)

- **Data:** Vizgen MERFISH Mouse Receptor Map (483 genes, 734k cells); flat CSV/HDF5
- **Config:** `data.vizgen_data_dir`
- **Deps:** same as Milestones 1–2 (`requirements.txt`; pandas CSV loader, no squidpy)
- **Goal:** Reproduce M1/M2 analyses; side-by-side comparison with Allen data for overlapping genes

### Milestone 4 — Zhuang (`04_zhuang_crossref.ipynb`)

- **Data:** `Zhuang-ABCA-1` … `Zhuang-ABCA-4` via same `AbcProjectCache` API
- **Goal:** Validate Allen MERFISH patterns; correlation scatter Allen vs Zhuang per region × cell type

See [`cursor_handover.md`](cursor_handover.md) for implementation sketches and Allen tutorial links.

---

## Known limitations and gotchas

1. **scRNA sub-regions are approximate** — CCF acronyms in config map to coarser dissection ROIs; warnings indicate when this happens.
2. **Large first-run downloads** — Plan disk space (~20 GB+ for cortex scRNA; ~50 GB if imputed MERFISH needed).
3. **MERFISH panel coverage** — Only ~500 genes measured; many receptors require imputed matrix or scRNA-only analysis.
4. **Per-gene MERFISH I/O** — Loading N genes opens the h5ad N times; slow for full 50-gene panel.
5. **Config vs notebook outputs** — Saved notebook outputs may reflect an older narrowed config; re-execute after YAML edits.
6. **Vizgen label transfer is approximate** — Vizgen cells get Allen `supertype` and CCF `brain_area` via kNN on overlapping genes, not native Vizgen annotations.
7. **figures/ in repo** — PNGs gitignored; `data/aggregated_scrna.parquet` may contain stale single-gene results until re-run.

---

## Suggested workflows

### Quick exploration (minimal download)

1. Narrow `receptor_query_config.yaml` to 1–2 genes and 1–2 broad regions (`STR`, `TH`)
2. Run M1 dev cell or `verify_setup.py` (without `--quick`)
3. Run M2 for genes in the 500-gene panel (check MERFISH gene list first)

### Full receptor survey

1. Use full config gene list; set `brain_areas` to regions of interest
2. Run `01_scrna_heatmaps.ipynb` — expect long runtime and large downloads for cortical regions
3. Run `02_merfish_spatial.ipynb` — set `use_imputed_merfish: false` first to avoid 50 GB download; enable for missing genes

### Cross-dataset validation (future)

1. Complete M1 + M2 on Allen data
2. Implement M3 (Vizgen) and M4 (Zhuang) per `cursor_handover.md`
3. Compare overlapping genes with correlation / side-by-side figures

---

## Document map

| Document | Audience | Content |
|----------|----------|---------|
| [README.md](README.md) | New users | Setup, notebook list, smoke test |
| **REPOSITORY_GUIDE.md** | Developers / analysts | Architecture, APIs, data flows, caveats |
| [cursor_handover.md](cursor_handover.md) | Implementation agents | Milestone specs, code sketches, Allen links |

---

*Generated from repository analysis. ABC Atlas manifest and download sizes reflect runs against `releases/20260415`.*
