# Methods — intrinsic physiology master build

## Executive summary

`build_intrinsic_master.py` (`data_parsing/intrinsic/`) transforms the Aug 2024 intrinsic-properties Excel workbook into structured outputs for downstream notebooks. The design goal is **faithful copy**: every value is parsed from source cells and attached to the correct `cell_id`; empty stays empty.

**Pipeline (Phases 1–4, all implemented)**

| Phase | When | What |
|-------|------|------|
| 1 — Build | Always | Parse 10 parameter sheets + metadata joins; deduplicate at 4 sig figs; write master + CSVs |
| 2 — Integrity | Pre-write | Re-read every parsed cell from Excel; assert pharmacology ⊆ control (cesium excepted) |
| 3 — Verification | Pre-write | Confirm cell-ID accounting, assembled values, label union, and full metadata joins |
| 4 — QC | Post-write | `Intrinsic_QC.xlsx`: stratified 3×IQR outliers, biological limits, log-scale flags |

Phases 2–3 **abort** on failure (no files written). Phase 4 is informational.

**Two analysis scopes**

- **`control_excitability`** — baseline measurements, one row per neuron (includes drug-sheet control blocks; cesium excluded).
- **`pharmacology_effect`** — effect blocks only, one row per `(cell_id, experiment)`.

**Key rules in one line each:** `excluded_in_May=1` cells dropped at parse; `exclude_flag` copied verbatim (informational only); duplicate headers disambiguated as `#1`/`#2`; overlapping cross-sheet values that disagree → `dup_conflict` + `duplicate_conflicts.csv`; numerics formatted to four significant figures on output.

**Default paths:** input `…/20260609_download_from_Drive/Intrinsic_properties_Analysis_Aug2024.xlsx`; output `…/physiology/restructured/`. Full provenance in `run_manifest.json`.

---

## Prime directive: copy-only, no recomputation

The build is a **faithful copy** task, not a compute task.

- **No hallucination.** Empty source cells stay empty in output.
- **No recomputation.** Values are not re-derived, unit-converted, or “cleaned” inside data cells.
- **Faithful parsing.** The engine uses `openpyxl` with `data_only=True`. Summary/stat columns (averages, SD, t-tests) are dropped at parse time so Excel error strings like `#DIV/0!` never propagate.
- **Exact ID ↔ data matching.** Each electrophysiology parameter is attached to the correct `cell_id` by construction of the sheet parsers; Phases 2–3 in the spec exist to prove this held.

Parameter columns are an intentional exception to byte-identical copy: numeric values are **coerced and formatted to four significant figures** for output (see [Parameter formatting](#parameter-formatting)). Deduplication also compares overlaps at four sig figs.

---

## Source workbook

**Default input:**  
`/Users/rancze/Documents/Data/expresso_data/physiology/20260609_download_from_Drive/Intrinsic_properties_Analysis_Aug2024.xlsx`

**Default output directory:**  
`/Users/rancze/Documents/Data/expresso_data/physiology/restructured/`

---

## Sheets included and excluded

### Parameter-data sheets (rows contribute measurements)

| Sheet | Routing |
|-------|---------|
| `V2M_L5`, `V1_L5`, `V1_L2-3`, `V2M_L2-3` | Control blocks → `control_excitability` |
| `All Analysed data` | Control blocks → `control_excitability` |
| `V2M_L5_2A_agonist` | Control → `control_excitability`; effect → `pharmacology_effect` (`5-HT2A agonist`) |
| `V2M_L5_1A_Antagonist` | Control → `control_excitability`; effect → `pharmacology_effect` (`5-HT1A antagonist`) |
| `V2M_L5_MDL` | Control → `control_excitability`; effect → `pharmacology_effect` (`MDL (5-HT2A antagonist)`) |
| `V2M_L5_TASK_Acidic_pH` | Control (pH 7.3) → `control_excitability`; effect (pH 6.1) → `pharmacology_effect` (`acidic pH (6.1)`) |
| `V2M_L5_Caesum` | Effect only → `pharmacology_effect` (`cesium (intracellular)`). **Excluded from `control_excitability`** — no control block. Sheet name **Caesum** is a student typo for cesium. |

### Metadata-only sheets (tags joined by cell ID; no parameter rows)

| Sheet | Output field |
|-------|----------------|
| `Assumed_tlx` | `assumed_type = Tlx` (mapped to `IT` on output) |
| `Assumed_PT_V2M` | `assumed_type = ET` |
| `cluster_analysis_res` | `physiological_cluster` (`IT` / `ET1` / `ET2`); may fill empty `assumed_type` |
| `All cells` | `exclude_flag`, `cre_label`, `axon`, `notes`, `time_from_5HT`, `excluded_in_May` filter |

The spec also lists `SC projecting cells` → `projection_target`; that sheet is **not** wired in the current build.

### Excluded entirely

Analysis, QC, summary, transposed, and duplicate sheets — e.g. `Excluded`, `Copy of *`, `ANR L5ET DATASET`, `ANR transposed`, `For_transposing`, `Transposed for plotting`, `ALL_CELLS_NEW`, `sig_diff_parameters_v1_vs_v2`, normality/ISSUE sheets, and all `cluster*` sheets except `cluster_analysis_res`. These never contribute parameter rows.

---

## Block detection and sheet-specific rules

### Standard area/layer sheets

- Cells start at the first `nm…` header column (~column I).
- **Condition rows** partition mega-blocks: `0` = control, `1` = effect; `2`, washout, and empty are ignored for control/effect routing.
- **Sub-sections** begin at column-A labels starting with `_`: `_IV`, `_short_depol`, `_EPSP`, `_crit_freq`, `_sag`, `_chirp`.
- The `_IV` header row itself carries data (`Vm_from_fit`) — it is not skipped.

### Drug sheets (`2A_agonist`, `1A_Antagonist`, `MDL`)

These introduce 2025 cells absent from area sheets. Each sheet has both a control block (`condition = 0`) and an effect block (`condition = 1`). Control feeds `control_excitability`; effect feeds `pharmacology_effect`.

### TASK (`V2M_L5_TASK_Acidic_pH`)

- Cell IDs from column B; duplicate headers are disambiguated (`#1` / `#2`).
- Control = section under **`pH 7.3`**; effect = section under **`pH 6.1`**. Section headers define control/effect — per-cell `0`/`1` condition rows are ignored (internally inconsistent).
- Only `_chirp`, `_sag`, `_EPSP` sub-sections are present. The pH 6.1 chirp block lacks a `_chirp` header and is inferred by matching parameter names against the pH 7.3 chirp block.
- Row 4 `exclude = 1` drops the cell from **both** output sheets (e.g. `nm2025_06_18_c1`).

### Cesium (`V2M_L5_Caesum`)

- Cells from column D (avg/SD columns B/C are dropped).
- No control block. The primary populated block is taken as the Cs⁺ effect.
- A sparse secondary block (~−56 mV, two cells) is skipped and noted as ambiguous.

### Parameter label union

All `(sub-section, label)` pairs across data sheets form the output parameter column set, ordered by sub-section (`_IV`, `_short_depol`, `_EPSP`, `_crit_freq`, `_sag`, `_chirp`) then alphabetically within section. Columns are named `section__label`.

Near-duplicate labels within the same sub-section are merged (string similarity ≥ 0.92, plus explicit merges in config). Every merge is logged in `run_manifest.json → label_merges`.

---

## Cell inclusion and exclusion

### `excluded_in_May`

Cells with `excluded_in_May = 1` on the **`All cells`** sheet are **dropped during parsing** and do not appear in any output file.

### TASK `exclude` flag

Cells with `exclude = 1` in row 4 of the TASK sheet are dropped from both scopes.

### `exclude_flag`

Copied verbatim from **`All cells`** column B, including non-numeric markers such as `?`. **Informational only** — not used to filter cells in downstream Expresso analysis.

---

## ID normalisation

- Strip whitespace; ensure `nm` prefix (`2025_06_19_c1` → `nm2025_06_19_c1`).
- **Duplicate headers:** when the same normalised ID appears in two columns on one sheet (known case: TASK `2025_06_19_c1`), disambiguate as `…_c1#1` / `…_c1#2` and record a `task_note`. These trigger **CRITICAL** terminal alerts and `duplicate_header_warnings` in the manifest.

---

## Control vs pharmacology routing (summary)

| Output | Scope | Contents |
|--------|-------|----------|
| `control_excitability` | Scope 1 | One row per unique neuron: metadata + union of control-block parameters. Includes drug-sheet control blocks. **Cesium excluded.** |
| `pharmacology_effect` | Scope 2 | One row per `(neuron × experiment)` with an effect block: same metadata and parameter columns plus `experiment`. |

Expected subset relationship (Phase 2): pharmacology cell IDs should ⊆ control IDs, **except** Cesium cells (effect-only).

---

## Metadata columns and column order

Metadata is joined by normalised `cell_id`. Column order in outputs:

**`control_excitability` and `pharmacology_effect`:**

1. `cell_id`
2. `region`
3. `areaCCF`
4. `layer`
5. `assumed_type`
6. `physiological_cluster`
7. `source_sheet`
8. `classic_burster`
9. `exclude_flag`
10. `cre_label`
11. `axon`
12. `notes` (merged from `All cells` note and comment, separated by ` | `)
13. `time_from_5HT`
14. `dup_conflict` (control sheet only; always empty/false on pharmacology rows)
15. `experiment` (pharmacology sheet only)
16. Parameter columns (`section__label`, ordered as above)

**`duplicate_conflicts`:**

`cell_id`, `region`, `areaCCF`, `layer`, `assumed_type`, `physiological_cluster`, `source_sheet`, `conflict_source_sheet`, then **only disagreeing parameter columns**.

### Provenance

| Field | Source |
|-------|--------|
| `source_sheet`, `layer` | Sheet name (`V1`/`V2M`, `L5`/`L2-3`) |
| `region` | Resolved broad atlas region (see below) |
| `areaCCF` | CCF subregion when resolved |
| `classic_burster` | Per-cell row on data sheet (`Classic burster? 0/1`) |
| `exclude_flag`, `cre_label`, `axon`, `time_from_5HT`, notes | `All cells` |
| `assumed_type` | `Assumed_tlx` / `Assumed_PT_V2M`; cluster fill; `PT` → `ET` |
| `physiological_cluster` | `cluster_analysis_res` column I (`IT`/`ET1`/`ET2`) |

### `region` vs `areaCCF`

There is **no separate `area` column** in output. Broad grouping is **`region`**: `VISp` or `V2M`. Source sheet labels `V1` are normalised to **`VISp`**.

**`areaCCF`** holds the finer CCF code when resolved: `VISp`, `VISam`, `VISpm`, `RSPagl`.

Resolution priority:

- **areaCCF:** morphology row on data sheet → `All cells` note → infer from sheet region (`VISp` when region is V1/VISp).
- **region:** derived from areaCCF (`VISp`→`VISp`, `VISam`/`VISpm`/`RSPagl`→`V2M`) → `All cells` area column → sheet region.

Builds report **`region_area_conflicts`** when sheet-derived region (`V1`/`V2M` from sheet name) disagrees with metadata/CCF-derived `region`. Sub-region vs broad-label pairs (e.g. VISpm with V2M) may appear but are often consistent; genuine mismatches (e.g. PTLp vs V2M) warrant manual review.

### `assumed_type`

- Source `Tlx` → output **`IT`**.
- `ET1` / `ET2` cluster labels map to **`ET`** for broad type; `physiological_cluster` keeps the verbatim `IT`/`ET1`/`ET2`.
- If `assumed_type` is empty but `physiological_cluster` is set, fill: `IT`→`IT`, `ET1`/`ET2`→`ET`.
- `PT` anywhere → `ET`.

---

## Parameter formatting

- Electrophysiology parameters are coerced to numeric where possible.
- Strings such as `"RMP cannot be determined"`, `#N/A`, `?`, `-`, and similar markers become **empty** in CSV/XLSX output.
- Numeric parameters are written with **four significant figures** (`format_param_for_output`).

This formatting applies to output files only; deduplication compares raw parsed values at four sig figs relative tolerance.

---

## Deduplication and `dup_conflict`

When the same `cell_id` appears on multiple included control sources, instances are grouped and overlapping parameters are compared at **four significant figures**.

- **All overlaps agree** → one canonical row. Canonical source is chosen by priority; non-overlapping parameters from lower-priority sheets are merged in.

**Dedup priority (highest first):**

`V2M_L5` → `V1_L5` → `V1_L2-3` → `V2M_L2-3` → `All Analysed data` → `V2M_L5_2A_agonist` → `V2M_L5_1A_Antagonist` → `V2M_L5_MDL` → `V2M_L5_TASK_Acidic_pH` → `V2M_L5_Caesum`

- **Any overlap disagrees** → canonical row kept in `control_excitability` with **`dup_conflict = TRUE`**; every source instance copied to `duplicate_conflicts`. Conflicting experimental values are a **data-quality red flag** for manual review.

Pharmacology rows are not deduplicated across conflicts in the same way; each `(cell_id, experiment)` takes the canonical instance by source priority if multiple exist.

---

## `duplicate_conflicts.csv` sparse format

The conflicts sheet lists **only parameter columns that actually disagree** between source instances — not the full parameter union. Each row includes metadata plus `conflict_source_sheet` identifying which workbook sheet the instance came from. Full per-cell detail (source pairs and raw disagreeing values) is also in `run_manifest.json → duplicate_conflicts`.

---

## `run_manifest.json`

Written beside the master workbook on every build. Key contents:

| Key | Meaning |
|-----|---------|
| `built_at` | Local timestamp |
| `git` | Repo root, commit hash (short and full), branch |
| `source_workbook` | Resolved path plus `created_at` / `modified_at` filesystem timestamps |
| `outputs` | Paths to XLSX and three CSVs |
| `counts` | Neurons, pharmacology rows, conflict cells/instances, `excluded_in_May` drops, parameter column count, cluster fills |
| `cluster_fill_column` | Column filled from cluster when `assumed_type` was empty (`assumed_type`) |
| `label_merges` | Parameter label normalisation log |
| `duplicate_header_warnings` | CRITICAL duplicate header list |
| `duplicate_conflicts` | Per-cell conflicting parameter keys and pair-wise detail |
| `region_area_conflicts` | Cells where sheet region ≠ metadata/CCF region |
| `phase2_verification` | Pre-write integrity gate results (byte-check count, cesium exceptions, pass/fail) |
| `phase3_verification` | Post-assembly verification (cell IDs, numerical/label/metadata fidelity) |
| `phase4_qc` | Outlier QC summary (flagged cell counts, QC workbook path) |
| `notes` | Human-readable explanations for warnings |

The CLI calls `print_run_alerts()` after the build: duplicate headers are shown first (loudest), then measurement conflicts, then region conflicts.

---

## Phase 2 — Pre-write integrity checks (implemented)

Before any CSV/XLSX is written, `intrinsic/verify.py` runs:

1. **Byte-equality re-read** — Every parsed `(sheet, row, col)` value is re-read from the workbook and compared with `==` to the in-memory copy. Count reported as `values_byte_checked` in `run_manifest.json → phase2_verification`.
2. **Pharmacology subset** — Every `pharmacology_effect` `cell_id` must appear in `control_excitability`, **except** cesium effect-only cells (listed in `cesium_effect_only_cells`).
3. **Accounting** — Orphan IDs on `ALL_CELLS_NEW` (informational; rescue not implemented), label merge count, conflict cells, cluster fills, `excluded_in_May` drops.

If any check fails, the build **aborts with exit code 1** and **no output files are written**. On success, terminal prints `Phase 2 integrity: PASSED`.

---

## Phase 3 — Post-build verification (implemented)

After parameter formatting (four sig figs) and before writing outputs, `intrinsic/verify_postbuild.py` verifies the assembled master rows:

1. **Cell-ID fidelity** — Every included source header ID appears in `control_excitability` or `pharmacology_effect`; no phantom IDs; TASK-dropped `nm2025_06_18_c1` absent; `#1`/`#2` disambiguated IDs present.
2. **Numerical fidelity** — Every parsed control/effect value matches the assembled row (after the same four-sig-fig formatting applied at write time). Duplicate-conflict non-canonical source instances are skipped.
3. **Label fidelity** — Output parameter column list equals the union of raw source labels after logged merges.
4. **Metadata fidelity** — Full join audit for `classic_burster`, `exclude_flag`, `cre_label`, `axon`, `time_from_5HT`, merged `notes`, and `physiological_cluster` against `All cells` / `cluster_analysis_res`.

If any check fails, the build **aborts** (no files written). On success, terminal prints `Phase 3 post-build verification: PASSED`. Results in `run_manifest.json → phase3_verification`.

---

## Phase 4 — Outlier QC workbook (implemented)

After the master files are written, `intrinsic/qc_outliers.py` produces **`Intrinsic_QC.xlsx`** (never merged into `Intrinsic_master.xlsx`). Three sheets:

| Sheet | Content |
|-------|---------|
| `flag_matrix` | Group mean/SD header rows, then all control and pharmacology rows; each parameter has a flag column and a twin `{parameter}__value` column (value shown only when flagged) |
| `param_summary` | Per group × parameter: n, mean, SD, median, Q1, Q3, min, max, IQR fences, outlier counts |
| `suspicious_cells` | Ranked list of cells/rows with the most flags |

**Stratification:** outliers are computed **within group only** — `VISp-L5`, `V2M-L5`, `VISp-L2/3`, `V2M-L2/3`. Pharmacology uses the same region–layer groups crossed with `experiment` (`region-layer|experiment`).

**Excluded from IQR pools:** rows with `exclude_flag = 1` remain in the master file but do not contribute to distribution statistics (`excluded_in_May` cells are already dropped at parse time).

**Flag matrix layout:** Before cell rows, two reference rows per group — `__GROUP_MEAN__` and `__GROUP_SD__` — populate the `{parameter}__value` twin columns with that group's mean and SD (flag columns empty). Cell rows place semicolon-separated flag codes in the parameter column and the outlying numeric value in the adjacent `__value` column.

**Flag types:**

| Code | Meaning |
|------|---------|
| `iqr` | Value outside the within-group 3×IQR fence (`Q1 − 3×IQR` or `Q3 + 3×IQR`) |
| `iqr+global` | Within-group 3×IQR outlier that is also outside the pooled global 3×IQR fence for that parameter (all groups in the same scope) |
| `log_iqr` | On log10-transformed values (Rin, rheobase, τ only; positive values), outside the within-group 3×IQR fence on the log scale |
| `log_iqr+global` | Log-scale 3×IQR outlier that is also outside the pooled global log-scale fence |
| `bio:Rin<0` | Input resistance below zero |
| `bio:tau<0` | Membrane time constant below zero |
| `bio:AP_peak<-20mV` | Action-potential peak more depolarized than −20 mV |
| `bio:Vm_rest>-30mV` | Resting membrane potential above −30 mV |
| `bio:rheobase_out_of_range` | Rheobase below 0 pA or above 2000 pA |
| `bio:AP_width<0.1ms` | AP width below 0.1 ms |
| `bio:sag_pct_outside_(0,100)` | Sag percentage not strictly between 0 and 100 |
| `bio:ratio_outside_(0,1)` | Adaptation or burst ratio not strictly between 0 and 1 |
| `bio:res_freq_out_of_range` | Resonance frequency below 0 Hz or above 200 Hz |

Semicolon-separated codes on one cell mean multiple rules fired (e.g. `iqr+global;bio:ratio_outside_(0,1)`). Canonical definitions are also written to `run_manifest.json → phase4_qc → flag_types` on each build.

Pharmacology flags use **absolute effect values**, not deltas. Summary in `run_manifest.json → phase4_qc`; terminal prints flagged cell/row counts.

---

## Cluster mapping (`cluster_analysis_res`)

- Columns from K onwards are disregarded.
- `cluster` (`IT`, `ET1`, `ET2`) is copied verbatim to `physiological_cluster`.
- Broad `assumed_type` is filled only when empty, using the mapping above.

---

## Cesium naming

The pharmacology sheet is named **`V2M_L5_Caesum`** in Excel (**Caesum** = student typo for cesium). Output experiment label: `cesium (intracellular)`. These cells appear only in `pharmacology_effect`.

---

## Related documents

- [README.md](README.md) — Quick start and output overview
- [Intrinsic_master_spec.md](Intrinsic_master_spec.md) — Full execution spec including Phases 2–4
- [METHODS_NOTES.md](METHODS_NOTES.md) — Original draft bullets (superseded for narrative by this file; kept for reference)
