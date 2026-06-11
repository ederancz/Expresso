# Intrinsic physiology data parsing

## Executive summary

This pipeline turns the Aug 2024 student Excel workbook (`Intrinsic_properties_Analysis_Aug2024.xlsx`) into analysis-ready tables for Expresso. It **copies** electrophysiology values and metadata as-is — no recomputation, imputation, or correction of experimental numbers.

**Run:** `conda activate expresso && python data_parsing/build_intrinsic_master.py`  
**Input → output:** physiology download folder → `…/physiology/restructured/`

**Main outputs**

| File | Role |
|------|------|
| `control_excitability.csv` | One row per neuron, baseline conditions (~97 cells) |
| `pharmacology_effect.csv` | One row per neuron × drug/manipulation |
| `duplicate_conflicts.csv` | Same `cell_id` on multiple sheets with disagreeing values |
| `Intrinsic_master.xlsx` | Excel version of the three sheets above |
| `Intrinsic_QC.xlsx` | Stratified outlier flags (informational; not in the master file) |
| `run_manifest.json` | Provenance, integrity results, conflict detail |

**Four phases:** (1) parse and assemble, (2) pre-write byte-equality check, (3) post-assembly verification, (4) QC outlier workbook. Phases 2–3 abort the build on failure; Phase 4 always runs after a successful write.

**Review after every build:** duplicate cell-ID headers (CRITICAL), measurement conflicts (`dup_conflict`), and region mismatches — all surfaced in the terminal and manifest.

---

## Prerequisites

- **Conda environment:** `expresso` (create with `conda env create -f environment.yml` from the repo root, then `conda activate expresso`).
- **openpyxl** — listed in `requirements.txt` (`openpyxl>=3.1`) and included in `environment.yml`.

## Default paths

Input and output live outside the repo under the Expresso data root:

| Role | Path |
|------|------|
| Source workbook | `/Users/rancze/Documents/Data/expresso_data/physiology/20260609_download_from_Drive/Intrinsic_properties_Analysis_Aug2024.xlsx` |
| Output directory | `/Users/rancze/Documents/Data/expresso_data/physiology/restructured/` |

Override either path with CLI flags (see below).

## How to run

From the repository root:

```bash
conda activate expresso
python data_parsing/build_intrinsic_master.py
```

Optional arguments:

```bash
python data_parsing/build_intrinsic_master.py \
  --source /path/to/Intrinsic_properties_Analysis_Aug2024.xlsx \
  --output-dir /path/to/restructured/
```

The script prints row counts and **Phase 2–4** pass/summary lines, then **terminal alerts** for data-quality issues (see below). If Phase 2 or Phase 3 fails, the build aborts and no output files are written. Full provenance is in `run_manifest.json`.

## Output files

| File | Description |
|------|-------------|
| `Intrinsic_master.xlsx` | Excel workbook with three sheets (header row frozen) |
| `control_excitability.csv` | CSV export of sheet 1 |
| `pharmacology_effect.csv` | CSV export of sheet 2 |
| `duplicate_conflicts.csv` | CSV export of sheet 3 |
| `Intrinsic_QC.xlsx` | Phase 4 QC workbook: flag matrix (with group mean/SD header rows), parameter summary, suspicious cells |
| `run_manifest.json` | Build metadata: git commit, counts, Phases 2–4 results, label merges, conflict detail |

## The three data sheets

**`control_excitability`** — One row per neuron under baseline (control) conditions. Includes area/layer sheets (`V1_L5`, `V2M_L2-3`, etc.), the pooled `All Analysed data` sheet, and **control blocks** from pharmacology sheets (drug experiments also record a pre-drug baseline). Cesium-only cells are excluded here because they have no control block.

**`pharmacology_effect`** — One row per `(neuron × experiment)` for drug or manipulation conditions: 5-HT2A agonist, 5-HT1A antagonist, MDL, acidic pH (TASK), and intracellular Cs⁺ (from the `V2M_L5_Caesum` sheet — “Caesum” is a typo in the source Excel). An `experiment` column identifies the condition.

**`duplicate_conflicts`** — When the same `cell_id` appears on multiple included sheets with **overlapping parameters that disagree beyond four significant figures**, the canonical row in `control_excitability` is flagged `dup_conflict = TRUE` and every conflicting source instance is copied here for manual review. Only parameter columns that actually disagree are populated (not the full parameter set).

## Terminal alerts

After a successful build, watch for banner messages:

- **CRITICAL — duplicate cell ID headers:** The same ID appears as two column headers on one sheet (almost always a student typo). Parsed as separate neurons with `#1` / `#2` suffixes. Review before analysis.
- **Duplicate measurement conflicts:** Same `cell_id` on multiple sheets with disagreeing values. Inspect `duplicate_conflicts.csv` and `run_manifest.json → duplicate_conflicts`.
- **Region conflicts:** Sheet-derived region disagrees with metadata/CCF-derived region for some cells (`run_manifest.json → region_area_conflicts`).

## Phase 4 QC — flag types

`Intrinsic_QC.xlsx` is informational only; it does not filter the master file. On **`flag_matrix`**, each parameter has a **twin column** (`parameter__value`) showing the numeric value when that cell is flagged. The sheet opens with **`__GROUP_MEAN__`** and **`__GROUP_SD__`** rows per stratified group (control: region–layer; pharmacology: region–layer|experiment) so flagged values can be read against group context.

Multiple flags on one value are semicolon-separated. Codes:

| Code | Meaning |
|------|---------|
| `iqr` | Outside within-group 3×IQR fence (Q1 − 3×IQR or Q3 + 3×IQR) |
| `iqr+global` | Within-group 3×IQR outlier **and** outside pooled global 3×IQR for that parameter |
| `log_iqr` | Log10-scale 3×IQR outlier (Rin, rheobase, τ; positive values only) |
| `log_iqr+global` | Log-scale outlier that is also outside the global log-scale fence |
| `bio:Rin<0` | Input resistance below zero |
| `bio:tau<0` | Membrane time constant below zero |
| `bio:AP_peak<-20mV` | AP peak more depolarized than −20 mV |
| `bio:Vm_rest>-30mV` | Resting Vm above −30 mV |
| `bio:rheobase_out_of_range` | Rheobase &lt; 0 or &gt; 2000 pA |
| `bio:AP_width<0.1ms` | AP width below 0.1 ms |
| `bio:sag_pct_outside_(0,100)` | Sag percentage not in (0, 100) |
| `bio:ratio_outside_(0,1)` | Adaptation or burst ratio not in (0, 1) |
| `bio:res_freq_out_of_range` | Resonance frequency &lt; 0 or &gt; 200 Hz |

Cells with `exclude_flag = 1` or `excluded_in_May = 1` on **`All cells`** are dropped at parse and do not appear in the master or QC output.

## Further reading

- **[METHODS.md](METHODS.md)** — Full methodology: sheet inclusion rules, ID normalisation, deduplication, metadata columns, and manifest format.
- **[Intrinsic_master_spec.md](Intrinsic_master_spec.md)** — Execution spec (Phases 1–4).
