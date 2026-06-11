# Methods notes — intrinsic physiology parsing

Draft bullets for README / METHODS (not user-facing yet).

## Cell inclusion

- Cells with `excluded_in_May = 1` in the source **All cells** sheet are **dropped during parsing** and do not appear in any output file.
- TASK sheet cells with `exclude = 1` (row 4) are also dropped.

## `exclude_flag`

- Copied verbatim from the source workbook (`All cells` column B), including non-numeric markers such as `?`.
- **Informational only** — not used to filter cells in downstream Expresso analysis.

## Duplicate measurements (`dup_conflict`)

- When the same `cell_id` appears on multiple included data sheets with **overlapping parameter values that disagree beyond 4 significant figures**, the build flags `dup_conflict = TRUE` on the canonical row in `control_excitability.csv` and writes all conflicting instances to `duplicate_conflicts.csv`.
- Full detail is recorded in `run_manifest.json` (git commit + per-cell conflict list).
- Conflicting experimental values are a **data-quality red flag** for manual review.

## Parameter formatting

- Electrophysiology parameter columns are coerced to numeric where possible; strings such as `"RMP cannot be determined"` become empty (NaN).
- Numeric parameters are written with **4 significant figures**.

## Cell typing

- Output column `assumed_type`: `Tlx` from source metadata is mapped to **`IT`**; `ET1`/`ET2` cluster labels map to **`ET`**.

## Area / region columns

- **`region`**: broad atlas grouping — **`VISp`** or **`V2M`**. Source sheet labels `V1` are normalised to `VISp`. This replaces the former separate `area` column.
- **`areaCCF`**: CCF subregion when resolved (`VISp`, `VISam`, `VISpm`, `RSPagl`).
- Builds report **`region_area_conflicts`** when sheet-derived region and metadata/CCF-derived region disagree.

## Duplicate conflicts file

- `duplicate_conflicts.csv` lists **only parameter columns that actually disagree** between source instances (not the full parameter set).
