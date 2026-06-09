# Expresso

Query and cross-validate gene expression (neurotransmitter **receptors** and intrinsic **excitability** ion channels) across mouse brain areas and cell types using the [Allen Brain Cell Atlas](https://alleninstitute.github.io/abc_atlas_access/) and complementary MERFISH datasets (Vizgen, Zhuang), within the Allen CCF + cell-type taxonomy framework.

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

For Milestone 3 (Vizgen), no extra pip packages are needed beyond `requirements.txt`. Download Vizgen CSVs and set `data.vizgen_data_dir` in the config.

Edit the config for genes and brain regions. The unified
[`query_config.yaml`](query_config.yaml) is the single source of truth: its
`gene_panel` is nested two levels — `category → family → genes` — with two
categories, `receptors` (neuromodulatory) and `excitability` (intrinsic ion
channels; see [`excitability_genes.md`](excitability_genes.md)). All notebooks
point at it via their `CONFIG_PATH` line.

To analyse only one category, comment out the other under `gene_panel`.
`load_config` also accepts a flat panel (`family → genes`) and the legacy
top-level keys `receptors`/`excitability`, so older configs still load unchanged.

## Notebooks

| Notebook | Description |
|----------|-------------|
| [`notebooks/01_scrna_heatmaps.ipynb`](notebooks/01_scrna_heatmaps.ipynb) | WMB-10Xv3 heatmaps by cell type × brain area |
| [`notebooks/02_merfish_spatial.ipynb`](notebooks/02_merfish_spatial.ipynb) | Allen MERFISH heatmaps (cell type × CCF brain area) |
| [`notebooks/03_vizgen_crossref.ipynb`](notebooks/03_vizgen_crossref.ipynb) | Vizgen MERFISH cross-ref vs Allen |
| [`notebooks/04_zhuang_crossref.ipynb`](notebooks/04_zhuang_crossref.ipynb) | Zhuang MERFISH cross-ref vs Allen |
| [`notebooks/05_synthesis.ipynb`](notebooks/05_synthesis.ipynb) | Cross-dataset evidence table, concordance & confidence tiers |

## Smoke test

```bash
python scripts/verify_setup.py --quick
```

Use `python scripts/verify_setup.py` (without `--quick`) to also test backed Drd2 loading from WMB-10Xv3-STR (large download).
