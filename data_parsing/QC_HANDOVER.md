# QC review handover

**To:** Hugh   
**From:** Expresso restructure pipeline  
**Purpose:** Review automated outlier flags on the intrinsic electrophysiology dataset before downstream analysis.

You know this dataset well. This note explains what I built, how to read the QC workbook, and what I need you to check.

---

## What you are receiving

| File | What it is |
|------|------------|
| **`Intrinsic_properties_Analysis_Aug2024.xlsx`** | Your original analysis workbook (unchanged source). Use this to look up raw values, notes, and sheet context. |
| **`Intrinsic_QC.xlsx`** | Automated QC output — flags suspicious values but **does not modify or delete anything** from the master data. |
| **`duplicate_conflicts.csv`** | **Cross-sheet duplicate disagreements** — same `cell_id` on multiple source sheets with conflicting measurements (see below). |

The restructured master tables were built by **copying** values from the source workbook (no re-analysis). **`control_excitability.csv`** holds one canonical baseline row per neuron; **`duplicate_conflicts.csv`** lists every conflicting source instance side-by-side for manual review.

---

## Why I am asking you to review

The pipeline flags values that are **statistically unusual within a group** or **biologically implausible by simple rules**. A flag is not a verdict — it means “look at this.” Many flags will be real biology (e.g. bursters, weak cells, drug effects); others may be entry errors, empty cells read as numbers, analysis errors, or ambiguous parameter labels.

We need your judgement on:

1. **QC flags** (`Intrinsic_QC.xlsx`) — keep, annotate, or exclude.
2. **Duplicate conflicts** (`duplicate_conflicts.csv`) — which sheet’s values are authoritative when the same cell appears in multiple places with **different numbers**.

---

## How to read `Intrinsic_QC.xlsx`

Three sheets:

### 1. `flag_matrix` (main review sheet)

**Top rows — group reference**  
Before any cell IDs, each stratified group has two rows:

- **`__GROUP_MEAN__`** — group average for each parameter (in the `…__value` columns)
- **`__GROUP_SD__`** — group standard deviation

Groups for **baseline (control)** data: `VISp-L5`, `V2M-L5`, `VISp-L2/3`, `V2M-L2/3`.  
For **pharmacology** rows, groups are `region-layer|experiment` (e.g. `V2M-L5|5-HT2A agonist`).

Use these rows to compare a flagged value against its peers in the same group.

**Cell rows**  
Each row is one neuron (control) or one neuron × drug condition (pharmacology). Columns:

- Metadata: `cell_id`, `scope`, `region`, `layer`, `group`, `experiment`, `total_flags`
- For each parameter: a **flag column** and a twin **`parameter__value`** column
  - Flag column: empty = no issue for that parameter; otherwise one or more codes (see below)
  - Value column: the numeric value **only when flagged** (empty if not flagged)

**Important:** Cells with `exclude_flag = 1` or `excluded_in_May = 1` on **`All cells`** are **omitted from the restructured master and QC files** entirely.

### 2. `param_summary`

Per group × parameter: sample size, mean, SD, quartiles, IQR fences, and counts of how many cells triggered each flag type. Useful if you want the distribution context without scrolling the full matrix.

### 3. `suspicious_cells`

Same flags as the matrix, but **ranked** — cells with the most flags first. Good starting point for review.

---

## What the flag codes mean

Flags can stack on one value, separated by `;` (e.g. `iqr+global;bio:ratio_outside_(0,1)`).

### Statistical flags (compare to group neighbours)

| Code | Plain meaning |
|------|----------------|
| **`iqr`** | Value falls outside the **within-group** “far tail” — below Q1 − 3×IQR or above Q3 + 3×IQR. We use 3×IQR (not the usual 1.5×) to avoid over-flagging skewed ephys data. |
| **`iqr+global`** | Same as `iqr`, **and** the value is also outside the pooled IQR fence across **all** groups in that scope (control or pharmacology). Stronger signal that the value is extreme relative to the whole dataset, not just its layer/area. |
| **`log_iqr`** | For **Rin, rheobase, and τ** only: the value is an outlier on a **log10 scale** within the group (these parameters are often right-skewed). |
| **`log_iqr+global`** | Log-scale outlier within group **and** in the global pooled log-scale distribution. |

**“Global”** here means: we first ask “is this odd within VISp-L5 (or whatever group)?” — if yes, we also ask “is it odd compared to every cell in the control (or pharmacology) pool for this parameter?” The `+global` tag marks values that fail both tests.

### Biological flags (`bio:…`)

Rule-of-thumb limits — not a full physiological model. They catch sign errors, wrong units, or values that are hard to interpret.

| Code | Rule |
|------|------|
| `bio:Rin<0` | Input resistance negative |
| `bio:tau<0` | Membrane time constant negative |
| `bio:AP_peak<-20mV` | AP peak more depolarized than −20 mV |
| `bio:Vm_rest>-30mV` | Resting Vm above −30 mV |
| `bio:rheobase_out_of_range` | Rheobase &lt; 0 or &gt; 2000 pA |
| `bio:AP_width<0.1ms` | AP width &lt; 0.1 ms |
| `bio:sag_pct_outside_(0,100)` | Sag percentage not between 0 and 100 |
| `bio:ratio_outside_(0,1)` | Adaptation or burst ratio not between 0 and 1 |
| `bio:res_freq_out_of_range` | Resonance frequency &lt; 0 or &gt; 200 Hz |

A **`bio:`** flag can appear alone (value fails the rule but may sit comfortably within the group distribution) or together with **`iqr`** (unusual *and* rule-violating).

---

## Duplicate cells and `duplicate_conflicts.csv`

Many neurons appear on **more than one sheet** (e.g. `V2M_L5` and **`All Analysed data`**, or an area sheet plus a drug sheet’s control block). That is expected. The pipeline **merges** them into one row in `control_excitability.csv` when values agree; when they **disagree beyond four significant figures**, it flags a **measurement conflict**.

### How merging works (dedup priority)

When the same `cell_id` has baseline data on multiple sheets, **overlapping parameters** are taken from the **highest-priority** sheet. Current order (highest first):

1. **`All Analysed data`**
2. `V2M_L5`
3. `V1_L5`
4. `V1_L2-3`
5. `V2M_L2-3`
6. Drug-sheet control blocks (`2A_agonist`, `1A_Antagonist`, `MDL`, `TASK`, …)

**Non-overlapping** parameters from lower-priority sheets are still merged in. **Region/layer** metadata come from the first area sheet that has them, even when parameters come from `All Analysed data`.

We configured **`All Analysed data` to win** on purpose — if you corrected values there after the area sheets were filled, the canonical row should now reflect the pooled sheet. **Conflicts in the CSV are exactly where that correction still doesn’t match the area sheet** (or another source).

### What lands in the main table vs the conflicts CSV

| Location | Meaning |
|----------|---------|
| **`control_excitability.csv`** | One **canonical** row per cell. Column **`dup_conflict = TRUE`** means at least one parameter disagreed across sources — the row still uses the priority rules above, not an average. |
| **`duplicate_conflicts.csv`** | **One row per conflicting source instance** for those cells. Only parameter columns that **actually disagree** are filled (sparse — not the full 65-parameter set). |

### Columns in `duplicate_conflicts.csv`

| Column | Meaning |
|--------|---------|
| `cell_id` | Neuron ID (same as source workbook, `nm` prefix) |
| `source_sheet` | Sheet of the **canonical** row in `control_excitability` (often `All Analysed data`) |
| **`conflict_source_sheet`** | Sheet **this row’s values** came from — compare rows sharing the same `cell_id` |
| `region`, `areaCCF`, `layer`, … | Metadata for that instance |
| Parameter columns | Populated **only** where this instance disagrees with another source |

**How to read it:** filter/sort by `cell_id`. You will see **2+ rows** per conflicted cell (e.g. one row with `conflict_source_sheet = V2M_L5`, another with `All Analysed data`), with the disagreeing parameters filled in each row. Open those same columns in the **source workbook** on both sheets and decide which value is correct.

### What we need from you on conflicts

For each conflicted `cell_id` (or at least the ones you recognise as your edits):

- **Confirm** which sheet should be authoritative (often `All Analysed data` if that’s where you consolidated fixes).
- **Fix the source workbook** so all sheets agree, **or** leave a clear note on `All cells` if one sheet is intentionally stale.
- After fixes, we rebuild; the cell should drop off `duplicate_conflicts.csv` and `dup_conflict` should become false.

**Latest build:** **23 cells** with conflicts, **46 instance rows** in the CSV (typically two sheets per cell).

**Not the same as:** duplicate **column headers** on TASK (`nm2025_06_19_c1#1` / `#2`) — that is a separate typo issue, not cross-sheet dedup.

---

## How to cross-check with your source workbook

1. Take `cell_id` from the QC sheet (same IDs as in your sheets, with `nm` prefix).
2. Open the matching **source sheet** (`source_sheet` is in the restructured CSVs if needed; pharmacology cells map to drug sheets or area sheets).
3. Confirm the flagged **parameter label** and raw cell value match what you expect.
4. Check **`All cells`** for `exclude_flag`, notes, and comments — you may have already documented why a cell is odd.

**Known source issues to keep in mind (already logged in the build):**

- **Duplicate column header** on TASK sheet: `2025_06_19_c1` appears twice → parsed as `nm2025_06_19_c1#1` and `#2`. Verify both columns are intentional.
- **`duplicate_conflicts.csv`:** cross-sheet measurement disagreements — review alongside QC (see section above).
- **`V2M_L5_Caesum`:** sheet name typo for cesium; cesium cells are effect-only (no baseline row in control table).
- **Excluded cells:** `exclude_flag = 1` or `excluded_in_May = 1` on **`All cells`** → dropped from restructured outputs (34 + 10 IDs in latest build).

---

## Suggested review workflow

**Duplicate conflicts (high priority if you edited `All Analysed data`):**

1. Open **`duplicate_conflicts.csv`** — sort by `cell_id`.
2. For each conflicted cell, compare rows differing in `conflict_source_sheet` (e.g. `V2M_L5` vs `All Analysed data`).
3. In the **source workbook**, check the listed parameter columns on both sheets.
4. Record which value is correct; update the workbook or note on `All cells`.

**QC outliers:**

1. Start with **`suspicious_cells`** in `Intrinsic_QC.xlsx` — work down from highest `total_flags`.
2. For each flagged cell, open **`flag_matrix`**, note the parameter and **`__value`**, and compare to **`__GROUP_MEAN__` / `__GROUP_SD__`** for that `group`.
3. In the **source workbook**, verify the underlying measurement and any sheet quirks (condition block, pH section, empty vs numeric).
4. Decide for each flag (or each cell):
   - **OK** — real biology or known protocol artefact; no action
   - **Annotation** — add or update note in `All cells` / lab records
   - **Exclude** — set `exclude_flag = 1` or `excluded_in_May = 1` on **`All cells`** (pipeline drops these on rebuild)

Please track decisions in whatever format works (marked-up CSV, comment list, or notes on `All cells`).

---

## What we are **not** asking you to do

- Re-run or rewrite the restructure pipeline
- Treat every flag as an error
- Change values in the source workbook without lab agreement (copy-only pipeline will not “fix” data on the next build)

---

## Quick numbers (latest build)

- ~**97** baseline (control) neurons, ~**36** pharmacology rows in master CSVs
- **23** cells with cross-sheet **`dup_conflict`**, **46** rows in **`duplicate_conflicts.csv`**
- ~**66** canonical baseline rows sourced from **`All Analysed data`** (dedup priority)
- ~**90+** control cells and ~**30** pharmacology rows with at least one QC flag
- **65** parameters tested; pharmacology QC uses **absolute effect values**, not deltas
- **10** IDs dropped via `excluded_in_May`, **34** via `exclude_flag = 1` on **`All cells`** (sets may overlap)

---

## Questions

If a flag code or group definition is unclear, or you find systematic false positives (e.g. one parameter name always flags bursters), note the pattern — we can tune rules or document exceptions for analysis.

Thank you for reviewing this; your domain knowledge is the step the automation cannot replace.
