# Expresso

Query receptor gene expression across mouse brain areas and cell types using the [Allen Brain Cell Atlas](https://alleninstitute.github.io/abc_atlas_access/).

## Python version

**Use Python 3.12** (3.11 also works). Python 3.13+ is not supported — many dependencies lack prebuilt wheels.

| File | Purpose |
|------|---------|
| [`environment.yml`](environment.yml) | Recommended: conda env with `python=3.12` |
| [`.python-version`](.python-version) | pyenv / local tooling |
| [`pyproject.toml`](pyproject.toml) | `requires-python = ">=3.11,<3.13"` |

## Setup (recommended)

```bash
conda env create -f environment.yml
conda activate expresso
```

To refresh an existing env:

```bash
conda env update -f environment.yml --prune
```

## Setup (pip only)

With Python 3.11 or 3.12 active:

```bash
pip install -r requirements.txt
```

For Milestone 3 (Vizgen), additionally:

```bash
pip install -r requirements-vizgen.txt
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
