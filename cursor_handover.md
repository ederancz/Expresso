# Handover: Mouse Brain Receptor Expression Notebook

## Goal

Build a Jupyter notebook that queries receptor gene expression across mouse brain
areas and cell types, using the Allen Brain Cell Atlas (ABC Atlas) as the primary
data source, with cross-reference to alternative datasets in later milestones.
Configuration (genes, brain regions, output settings) is read entirely from
`receptor_query_config.yaml` — no hard-coded gene names or region lists anywhere
in the notebook code.

---

## Repository layout

```
project/
├── receptor_query_config.yaml   # runtime config — edit to change genes/regions
├── notebooks/
│   ├── 01_scrna_heatmaps.ipynb  # Milestone 1
│   ├── 02_merfish_spatial.ipynb # Milestone 2
│   ├── 03_vizgen_crossref.ipynb # Milestone 3 (future)
│   └── 04_zhuang_crossref.ipynb # Milestone 4 (future)
├── src/
│   ├── config.py                # load + validate YAML config
│   ├── data_loaders.py          # ABC Atlas and Vizgen I/O helpers
│   ├── plotting.py              # heatmap + spatial plot functions
│   └── utils.py                 # gene-panel intersection, normalisation helpers
└── figures/                     # auto-created by notebooks
```

---

## Data sources

### Milestone 1 & 2 — Allen Brain Cell Atlas (primary)

- **Python package:** `abc_atlas_access`
  ```
  pip install git+https://github.com/alleninstitute/abc_atlas_access.git
  ```
- **Data location:** AWS S3 public bucket `arn:aws:s3:::allen-brain-cell-atlas`
  (no account needed; `AbcProjectCache` handles download/caching)
- **Key datasets:**

| ID | Type | Cells | Genes | Notes |
|---|---|---|---|---|
| `WMB-10Xv3` | scRNA-seq | ~4M | ~32k | Whole mouse brain |
| `MERFISH-C57BL6J-638850` | MERFISH | ~4.3M | ~500 | Spatial + CCF coords |
| `MERFISH-C57BL6J-638850-imputed` | MERFISH | ~4.3M | ~8k | Imputed from scRNA ref |
| `Zhuang-ABCA-1` through `-4` | MERFISH | ~4 replicates | ~1,100 | Harvard lab, independent |

- **Cell type hierarchy:** class → subclass → supertype → cluster
  (34 / 338 / 1,201 / 5,322 groups)
- **Brain area annotation:** CCF v3 parcellation attached to every cell

### Milestone 3 — Vizgen MERFISH Mouse Receptor Map (cross-reference)

- 483-gene panel (GPCRs, RTKs, canonical markers); 734k cells; 3 coronal slices × 3 replicates
- Flat files (CSV + HDF5), no special API; download from: https://info.vizgen.com/mouse-brain-data
- Tutorial available in squidpy: https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_vizgen.html

### Milestone 4 — Zhuang MERFISH (cross-reference)

- Already accessible via `AbcProjectCache` as `Zhuang-ABCA-{1,2,3,4}`
- ~1,100-gene panel; 4 whole-brain replicates

---

## Milestone 1 — scRNA-seq heatmaps (`01_scrna_heatmaps.ipynb`)

### Objective
For each receptor in config, compute mean log2(CPM+1) expression per
`{cell_type_level} × brain_area` combination. Produce:
- One heatmap per receptor family (rows = cell types, columns = brain areas)
- One combined heatmap (all receptors × top-N cell types, clustered)

### Key implementation notes

**Partial loading pattern** — do NOT load the full expression matrix:
```python
# AbcProjectCache splits expression by anatomical package; load only needed genes
import anndata as ad

def load_genes_for_regions(cache, dataset_id, genes, brain_areas, cell_meta):
    """Load a genes x cells slice without reading the full matrix."""
    # filter cell metadata to target regions first
    mask = cell_meta['region_of_interest_acronym'].isin(brain_areas)
    target_cells = cell_meta[mask].index

    # expression files are split by anatomical package; iterate only relevant ones
    frames = []
    for pkg in cache.get_expression_matrix_packages(dataset_id):
        ad_file = cache.get_expression_matrix_path(dataset_id, pkg)
        # backed='r' avoids loading the full matrix into RAM
        adata = ad.read_h5ad(ad_file, backed='r')
        cells_in_pkg = adata.obs_names.intersection(target_cells)
        if len(cells_in_pkg) == 0:
            continue
        genes_in_pkg = [g for g in genes if g in adata.var_names]
        if len(genes_in_pkg) == 0:
            continue
        subset = adata[cells_in_pkg, genes_in_pkg].to_memory()
        frames.append(subset)
    return ad.concat(frames) if frames else None
```

**Aggregation:**
```python
import pandas as pd, numpy as np

def aggregate_expression(adata, cell_meta, cell_type_col, region_col):
    df = pd.DataFrame(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                      index=adata.obs_names, columns=adata.var_names)
    df[cell_type_col] = cell_meta.loc[df.index, cell_type_col]
    df[region_col]    = cell_meta.loc[df.index, region_col]
    return df.groupby([cell_type_col, region_col]).mean()
```

**Gene panel intersection:** Always check which requested genes are actually present
in the dataset; warn (don't crash) for missing ones. For scRNA-seq, all ~32k genes
are available.

### Expected outputs
- `figures/heatmap_{family}.png` per receptor family
- `figures/heatmap_combined.png`
- `data/aggregated_scrna.parquet` (if `save_processed_data: true`)

---

## Milestone 2 — MERFISH spatial maps (`02_merfish_spatial.ipynb`)

### Objective
For each receptor, plot spatial expression across the whole brain (x/y/z CCF
coordinates), coloured by expression level. Produce:
- Per-receptor scatter plots (coronal, sagittal projections)
- Overlay with CCF parcellation boundaries (optional, if shapefile available)
- Side-by-side: raw MERFISH panel genes vs imputed genes (where applicable)

### Key implementation notes

**Gene availability:** The MERFISH panel has ~500 genes. For genes outside the panel,
load from the imputed matrix (`MERFISH-C57BL6J-638850-imputed`) if
`use_imputed_merfish: true` in config. Flag imputed genes clearly in plot titles.

**Spatial loading pattern:**
```python
# Cell metadata (with x, y, z, CCF annotation) is a lightweight CSV — load fully
cell_meta = cache.get_metadata_dataframe('MERFISH-C57BL6J-638850', 'cell_metadata')

# Load expression for one gene at a time to minimise RAM
def load_single_gene_merfish(cache, gene, cell_meta, dataset_id):
    panel_genes = cache.get_gene_list(dataset_id)
    if gene in panel_genes:
        src = dataset_id
    elif config['data']['use_imputed_merfish']:
        src = dataset_id + '-imputed'
        # check gene exists in imputed panel
    else:
        return None, 'not_in_panel'
    adata = ad.read_h5ad(cache.get_expression_matrix_path(src), backed='r')
    expr = np.asarray(adata[:, gene].X).ravel()
    return pd.Series(expr, index=adata.obs_names), src
```

**Spatial plot function:**
```python
def plot_spatial(coords_df, expression, gene, projection='coronal',
                 cmap='magma', source_label=''):
    """
    coords_df: DataFrame with columns x, y, z (CCF µm)
    projection: 'coronal' (x,y), 'sagittal' (z,y), 'axial' (x,z)
    """
    axes = {'coronal': ('x','y'), 'sagittal': ('z','y'), 'axial': ('x','z')}
    ax1, ax2 = axes[projection]
    fig, ax = plt.subplots(figsize=config['output']['figsize_spatial'])
    sc = ax.scatter(coords_df[ax1], coords_df[ax2],
                    c=expression, cmap=cmap, s=0.3, alpha=0.6,
                    vmin=0, vmax=np.percentile(expression[expression>0], 99))
    plt.colorbar(sc, ax=ax, label='log2(CPM+1)')
    title = f'{gene} — MERFISH {projection}'
    if source_label == 'imputed':
        title += ' [imputed]'
    ax.set_title(title)
    ax.set_aspect('equal')
    ax.invert_yaxis()
```

### Expected outputs
- `figures/spatial_{gene}_{projection}.png` per gene × projection
- `figures/spatial_panel_{family}.png` (multi-panel grid per receptor family)

---

## Milestone 3 — Vizgen cross-reference (`03_vizgen_crossref.ipynb`) *[future]*

### Objective
Load the Vizgen MERFISH Mouse Receptor Map (483-gene panel), reproduce the
heatmap and spatial analyses from Milestones 1–2, then produce a side-by-side
comparison figure for genes that overlap between datasets.

### Implementation sketch
- Data location: `config['data']['vizgen_data_dir']` (user downloads flat files)
- Use `squidpy` or plain pandas/numpy to parse Vizgen cell-by-gene matrix
- Map Vizgen spatial coordinates to CCF if possible (affine transform provided
  in some publications); otherwise use native Vizgen coordinates
- Gene overlap: `set(config_genes) ∩ set(vizgen_panel)` — currently 483 genes

---

## Milestone 4 — Zhuang MERFISH cross-reference (`04_zhuang_crossref.ipynb`) *[future]*

### Objective
Access the 4 Zhuang whole-brain MERFISH replicates (`Zhuang-ABCA-{1,2,3,4}`)
via `AbcProjectCache` and validate expression patterns observed in Milestone 2.

### Implementation sketch
- Same `AbcProjectCache` API as Milestone 2; dataset IDs: `Zhuang-ABCA-1` … `Zhuang-ABCA-4`
- Tutorial to reference: https://alleninstitute.github.io/abc_atlas_access/notebooks/zhuang_merfish_tutorial.html
- ~1,100-gene panel; check overlap with config gene list
- Produce correlation scatter (Allen MERFISH vs Zhuang MERFISH) per region × cell type

---

## Shared utilities (`src/`)

### `config.py`
```python
import yaml
from pathlib import Path

def load_config(path='receptor_query_config.yaml'):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    # Flatten receptor families into a single list with family labels
    genes = {}
    for family, glist in cfg['receptors'].items():
        for g in glist:
            genes[g] = family
    cfg['_genes_flat'] = genes          # {gene: family}
    cfg['_all_genes']  = list(genes)    # ordered list
    return cfg
```

### `data_loaders.py`
- `get_abc_cache(config)` — initialise `AbcProjectCache` from `cache_dir`
- `load_cell_metadata(cache, dataset_id, brain_areas)` — returns filtered DataFrame
- `load_expression_subset(cache, dataset_id, genes, cell_ids)` — backed h5ad partial load
- `check_gene_availability(cache, dataset_id, genes)` — returns `{gene: 'present'|'imputed'|'missing'}`

### `plotting.py`
- `plot_heatmap(agg_df, genes, title, config)` — seaborn clustermap wrapper
- `plot_spatial(cell_meta, expression, gene, projection, config)`
- `plot_family_panel(results, family, config)` — multi-receptor figure grid

---

## Dependencies

```
# requirements.txt
git+https://github.com/alleninstitute/abc_atlas_access.git
anndata>=0.10
pandas>=2.0
numpy>=1.26
scipy>=1.12
matplotlib>=3.8
seaborn>=0.13
squidpy>=1.4          # for Vizgen milestone
pyyaml>=6.0
pyarrow>=14.0         # for parquet cache
tqdm
```

---

## Reference notebooks (read before implementing)

Study these Allen Institute notebooks before writing any I/O code — they show
the exact `AbcProjectCache` API and file structure:

1. Getting started: https://alleninstitute.github.io/abc_atlas_access/notebooks/getting_started.html
2. 10x gene expression (part 1): https://alleninstitute.github.io/abc_atlas_access/notebooks/10x_snRNASeq_tutorial_part_1.html
3. MERFISH spatial (part 1): https://alleninstitute.github.io/abc_atlas_access/notebooks/merfish_tutorial_part_1.html
4. MERFISH imputed genes: https://alleninstitute.github.io/abc_atlas_access/notebooks/merfish_imputed_genes_example.html
5. CCF registration: https://alleninstitute.github.io/abc_atlas_access/notebooks/merfish_ccf_registration_tutorial.html


---


