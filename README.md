# Expresso

Query receptor gene expression across mouse brain areas and cell types using the [Allen Brain Cell Atlas](https://alleninstitute.github.io/abc_atlas_access/).

## Setup

```bash
pip install -r requirements.txt
```

Edit [`receptor_query_config.yaml`](receptor_query_config.yaml) for genes and brain regions.

## Notebooks

| Notebook | Description |
|----------|-------------|
| [`notebooks/01_scrna_heatmaps.ipynb`](notebooks/01_scrna_heatmaps.ipynb) | WMB-10Xv3 heatmaps by cell type × brain area |
| [`notebooks/02_merfish_spatial.ipynb`](notebooks/02_merfish_spatial.ipynb) | MERFISH spatial maps (CCF coordinates) |
| `notebooks/03_vizgen_crossref.ipynb` | TODO — Milestone 3 |
| `notebooks/04_zhuang_crossref.ipynb` | TODO — Milestone 4 |

## Smoke test

```bash
python scripts/verify_setup.py --quick
```

Use `python scripts/verify_setup.py` (without `--quick`) to also test backed Drd2 loading from WMB-10Xv3-STR (large download).
