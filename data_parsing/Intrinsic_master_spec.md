# Execution spec — Intrinsic-properties master build + QC

**Handover document for implementation (Cursor).**

---

## 0. Prime directive — data fidelity (non-negotiable)

This is a pure **copy** task, not a compute task. The single most important
requirement:

- **No hallucination.** Never invent, infer, interpolate, or "correct" a value.
  If a source cell is empty, the output is empty.
- **Faithful numerical parsing.** Every numeric value is carried byte-for-byte
  from source to output. Never recompute, round, reformat, or change type
  (a string stays a string; a number stays a number).
- **Never modify, only copy.** No unit conversion, no normalisation of values,
  no cleanup of stray characters inside data cells.
- **Perfect ID ↔ data matching.** The join between a cell ID and its
  parameters / comments / metadata must be exact. A misaligned column or an
  off-by-one row that attaches the wrong value to the wrong cell is the worst
  possible failure mode and must be impossible by construction (Phase 2 + 3
  exist to prove this).

Engine: `openpyxl`, `data_only=True`. Summary columns (avg/SD/ttest) are dropped,
so error strings like `#DIV/0!` never propagate.

Output: workbook `Intrinsic_master.xlsx` with 3 sheets, plus a CSV export of each
sheet. Freeze the header row. Rows = cells. A **separate** QC workbook is produced
in Phase 4 — do **not** add QC sheets to `Intrinsic_master.xlsx`.

---

## Phase 1 — Build

### Sheet inclusion

**Parameter-data sources:**
`V2M_L5`, `V1_L5`, `V1_L2-3`, `V2M_L2-3`, `All Analysed data`,
`V2M_L5_2A_agonist`, `V2M_L5_1A_Antagonist`, `V2M_L5_MDL`,
`V2M_L5_TASK_Acidic_pH`, `V2M_L5_Caesum`.

**Metadata-only sources** (no parameter rows; supply tags + feed identity check):
- `Assumed_tlx` → `assumed_type = Tlx`
- `Assumed_PT_V2M` → `assumed_type = ET`
- `SC projecting cells` → `projection_target = SC`
- `cluster_analysis_res` → `physiological_cluster` (see **Cluster mapping** below)

**Excluded entirely:**
`Excluded`, all `Copy of *`, `ANR L5ET DATASET`, `ANR transposed`,
`For_transposing`, `Transposed for plotting`, `ALL_CELLS_NEW`,
`sig_diff_parameters_v1_vs_v2`, `compare <5min to 9-10 in 5HT V1`,
all QC / summary / normality / ISSUE sheets, all `cluster*` sheets
**except `cluster_analysis_res`**, and any other analysis / comparison /
transposed sheet. (`ANR` and `ALL_CELLS_NEW` may be cross-checked for identity
but never contribute parameter rows.)

**Orphan-ID rescue:** if any excluded/transposed sheet (e.g. `ALL_CELLS_NEW`)
contains a cell ID that appears on **no** included data sheet and is not
explicitly marked excluded anywhere, include that cell and report it. (In the
known dataset `ALL_CELLS_NEW` had zero unique IDs — expect none, but verify.)

### Cluster mapping (`cluster_analysis_res`)

- **Disregard columns from K onwards** on this sheet.
- The `cluster` column separates one IT type and two physiological ET types:
  `IT`, `ET1`, `ET2`.
- Copy the `cluster` value **verbatim** into a new output metadata column
  **`physiological_cluster`** (`IT` / `ET1` / `ET2`), joined by normalised ID.
  This column is carried into **both** output sheets' CSVs.
- For cells listed here, if the broad cell-type field is **unpopulated**, fill it
  using: `IT → IT`, `ET1 → ET`, `ET2 → ET`. The verbatim `physiological_cluster`
  is never collapsed — only the broad type is filled.
  - ⚠️ **CONFIRM which field is "the cell-type column."** Best current
    interpretation: the broad type is `assumed_type` (currently `Tlx`/`ET`,
    with `PT → ET` already normalised). If the dataset has a distinct explicit
    `cell_type`/`IT-ET` column, fill that one instead. Do not guess silently —
    surface which column was filled in the build report.

### Cesium handling (effect-only)

`V2M_L5_Caesum` **is included in `pharmacology_effect` for completeness**, as the
`Cs⁺ (intracellular)` experiment. Cesium cells are **excluded from
`control_excitability`** (they have no control block). This restores the original
effect-only routing; the subset assertion (Phase 2) lists Cesium as the expected
exception.

### ID normalisation

Strip whitespace; ensure `nm` prefix (odd sheets use `2025_..._cN` → prepend
`nm`). TASK has a duplicated header `2025_06_19_c1` (two columns) — keep both,
disambiguate as `..._c1#1` / `..._c1#2`, and flag in a note column (`task_note`).

### Block detection

**Standard sheets:** partition by condition rows. Mega-block with per-cell
condition `= 0` is **control**; `= 1` is **effect**; `2`/`Washout`/empty ignored.
Sub-sections begin at col-A labels starting `_` (`_IV`, `_short_depol`, `_EPSP`,
`_crit_freq`, `_sag`, `_chirp`). Map `(sub-section, param-label) → per-cell value`.
Cells start at the first `nm…` header column (≈ col I).
**Note:** the `_IV` header row itself carries data (`Vm_from_fit`) — do not skip it.

**Drug sheets** (`2A_agonist`, `1A_Antagonist`, `MDL`, `TASK`): these contain
entirely new 2025 cells absent from the area / All-Analysed sheets. Each has
**both** a control (condition=0) block and an effect block. The control block
**feeds `control_excitability`** (otherwise every drug cell violates the
scope-2 ⊆ scope-1 assertion). The effect block feeds `pharmacology_effect`.

**TASK (special-case):** cells from col B. Control = section under the `pH 7.3`
header; effect = section under the `pH 6.1` header. Use the **section headers** as
the control/effect signal; **ignore** the per-cell 0/1 condition rows (internally
inconsistent). Only `_chirp` / `_sag` / `_EPSP` present. The `pH 6.1` chirp block
has **no `_chirp` sub-section header** — infer it as `_chirp` by parameter-name
matching against the `pH 7.3` chirp block. TASK row-4 exclusion flag: any cell with
`exclude = 1` (e.g. `nm2025_06_18_c1`) is **dropped entirely** from both scopes.

**Cesium (special-case):** cells from col D (avg/SD in B/C, dropped). No control.
Take the primary populated block as the `Cs⁺` effect → scope 2 only. Flag/skip the
sparse secondary (−56 mV, 2-cell) block as ambiguous in a note column.

### Parameter union

Ordered union of `(sub-section, canonical-label)` across all data sheets.
Auto-merge near-identical labels **within the same sub-section** (e.g.
`First AP peak in RecordA )mV)1` → `First AP peak in RecordA (mV)`;
`…RecordA (mV)1` variants). Build the merge map by within-section string
similarity; **log every merge** for review. Missing params → blank.

### Metadata columns (scope 1, carried to scope 2), joined by normalised ID

- From sheet name: `source_sheet`, `region` (V1/V2M), `layer` (L2-3/L5).
- From sheet per-cell rows: `classic_burster` (Classic burster? row),
  `area_morph` (morphology/area row, e.g. VISpm).
- From `All cells` (keyed by ID): `exclude_flag`, `area_meta`, `layer_meta`,
  `cre_label` (Colgalt+), `axon`, `note`, `comment`, `excluded_in_May`,
  `time_from_5HT`.
- From metadata-only sheets: `assumed_type` (Tlx/ET), `projection_target` (SC),
  `physiological_cluster` (IT/ET1/ET2).
- Keep `area_morph` and `area_meta` as **separate** columns; add
  `area_mismatch = TRUE` where both present and differ.
  **Note:** sub-region vs broader-label pairs (VISam/VISpm vs V2M) are
  *consistent*, not genuine mismatches — they will be flagged but are expected;
  genuine mismatches (e.g. PTLp vs V2M) are also flagged for manual review.
- Normalise any `PT → ET` wherever it appears.

### Outputs

- **`control_excitability`** — one row per unique neuron (max coverage), metadata
  columns + control-block parameter union. Includes drug-sheet control blocks.
  **Cesium excluded.**
- **`pharmacology_effect`** — one row per `(neuron × experiment)` with a
  condition-1/effect block, same parameter columns + an `experiment` column:
  - `*_2A_agonist` → `5-HT2A agonist`
  - `*_1A_Antagonist` → `5-HT1A antagonist`
  - `*_MDL` → `MDL (5-HT2A antagonist)`
  - `*_TASK_Acidic_pH` → `acidic pH (6.1)`
  - `*_Caesum` → `Cs⁺ (intracellular)`
- **`duplicate_conflicts`** — for cells whose overlapping parameter values disagree
  beyond 4 significant figures across sources: copy all instances as separate rows,
  tagged with `source_sheet`. Also keep the canonical instance in sheet 1/2 with
  `dup_conflict = TRUE`.

### Dedup rule

Group data-sheet instances by normalised ID. Compare overlapping params at 4 sig
figs (relative tolerance, to absorb rounding). All agree → keep canonical by
priority `V2M_L5 / V1_L5 / V1_L2-3 / V2M_L2-3 > All Analysed data > drug sheets`;
merge metadata. Any disagree → route to `duplicate_conflicts`.

---

## Phase 2 — Integrity checks (before write)

- **Re-read each copied value from source and assert byte-equality** — guarantee
  no mutation.
- Assert `pharmacology_effect` IDs ⊆ `control_excitability` IDs. **Expected
  exception: Cesium cells** (effect-only, no control) — list them; any *other*
  violation is an error.
- Report counts: neurons in scope 1; `(neuron × experiment)` rows in scope 2;
  conflict cells; label merges performed; orphan IDs rescued; cluster field-fills
  applied (and which column was filled).

---

## Phase 3 — Post-build verification (full, not spot-check)

This phase exists to *prove* the Prime Directive held.

- **Cell-ID fidelity:** every cell ID from every source row-1 header appears in
  output (no silent drops/additions); `nm`-prefix normalisation consistent;
  TASK `#1`/`#2` correct; dropped TASK cell (`nm2025_06_18_c1`) absent.
- **Numerical fidelity:** for each source data sheet × canonical section, sample
  cells and compare raw values to master; flag any mismatch beyond
  floating-point identity. Pay attention to the `_IV` header row (`Vm_from_fit`)
  and to any `None` in output (confirm source was truly empty vs a dropped value).
- **Label fidelity:** output header list == union of raw source labels after the
  logged merges; no label silently collapsed or split.
- **Metadata fidelity (FULL join, not sample):** cross-check `classic_burster`,
  `area_morph` against source sheet rows; cross-check `area_meta`, `layer_meta`,
  `cre_label`, `axon`, `time_from_5HT` against `All cells` for **every** cell;
  cross-check `physiological_cluster` against `cluster_analysis_res` for every
  cluster-listed cell.

---

## Phase 4 — Outlier analysis (separate QC workbook, e.g. `Intrinsic_QC.xlsx`)

- **Stratification:** flag outliers **within group only** — `V1-L5`, `V2M-L5`,
  `V1-L2/3`, `V2M-L2/3` computed separately. Do **not** pool across layers/regions.
  Note where a within-group outlier is also a global outlier.
- **Excluded cells:** cells with `exclude_flag = 1` **or** `excluded_in_May = 1`
  are **excluded from the outlier distributions** (they remain in the master file
  but must not skew the IQRs).
- **Pharmacology:** report outliers on **absolute values** in the drug condition
  (not deltas).

For each parameter column, per group:
1. Basic stats: `n, mean, SD, median, Q1, Q3, min, max`.
2. Statistical outliers: outside `Q1 − 3×IQR` or `Q3 + 3×IQR` (conservative Tukey
   fence; biological data are right-skewed, so 1.5× over-flags).
3. Biological plausibility hard limits (derive the full list from parameter names;
   examples: `Rin<0`, `τ<0`, `AP peak < −20 mV`, `Vm_rest > −30 mV`,
   `rheobase < 0 or > 2000 pA`, `AP width < 0.1 ms`, `sag ratio outside (0,1)`,
   `resonance freq < 0 or > 200 Hz`).
4. For skewed parameters (`Rin`, `rheobase`, `τ`): additionally test log-normality
   and flag outliers on the log scale.

Produce: a `cell × parameter` flag matrix, a per-parameter summary (per group),
and a ranked "most suspicious" cell list.

---

## Open item for the implementer to confirm

1. **Cluster cell-type fill target** — which existing column is "the cell-type
   column" to fill from `IT/ET1→IT/ET` (best guess: `assumed_type`). See the
   Cluster mapping section. Surface the chosen column in the build report.
