#!/usr/bin/env python3
"""One-off schema audit for control_excitability.csv (M6 phase 0)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

DEFAULT_PATH = Path(
    "/Users/rancze/Documents/Data/expresso_data/physiology/restructured/control_excitability.csv"
)


def main(path: Path) -> None:
    df = pd.read_csv(path, low_memory=False)
    print("=" * 60)
    print("FILE:", path)
    print("SHAPE:", df.shape)
    print("=" * 60)

    meta_cols = [c for c in df.columns if not c.startswith("_")]
    print(f"\nMETADATA COLUMNS ({len(meta_cols)}):")
    for c in meta_cols:
        nnull = df[c].isna().sum() + (df[c].astype(str).str.strip() == "").sum()
        print(f"  {c!r}: {df[c].nunique()} unique, ~{nnull} empty")

    chirp_cols = [c for c in df.columns if c.startswith("_chirp")]
    print(f"\nCHIRP COLUMNS ({len(chirp_cols)}):")
    for c in chirp_cols:
        print(f"  {c}")

    for col in ["region", "layer", "assumed_type", "classic_burster", "exclude_flag",
                "excluded_in_May", "area_mismatch", "dup_conflict"]:
        if col in df.columns:
            print(f"\n--- {col} ---")
            print(df[col].value_counts(dropna=False).to_string())

    print("\n--- area_morph (top 15) ---")
    print(df["area_morph"].value_counts(dropna=False).head(15).to_string())

    print("\n--- area_meta ---")
    print(df["area_meta"].value_counts(dropna=False).to_string())

    print("\n--- source_sheet (top 10) ---")
    print(df["source_sheet"].value_counts(dropna=False).head(10).to_string())

    rf = "_chirp__Res. freq. (Hz)"
    ri = "_chirp__Res. imp. mag. (MOhm)"
    for c in [rf, ri]:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            print(f"\n--- {c} ---")
            print(f"  non-null: {s.notna().sum()} / {len(df)}")
            if s.notna().any():
                print(s.describe().to_string())

    print("\n--- region × layer ---")
    print(pd.crosstab(df["region"], df["layer"], dropna=False).to_string())

    sub = df[df["assumed_type"].notna() & (df["assumed_type"].astype(str).str.strip() != "")]
    print("\n--- region × assumed_type (typed cells only, n=%d) ---" % len(sub))
    print(pd.crosstab(sub["region"], sub["assumed_type"], dropna=False).to_string())

    # Proposed brain_area_group mapping
    def map_group(r: str) -> str:
        r = str(r).strip()
        if r == "V1":
            return "VISp"
        if r == "V2M":
            return "V2M"
        return r or "UNKNOWN"

    df["_brain_area_group"] = df["region"].map(map_group)

    def map_coarse(row) -> str:
        layer = str(row.get("layer", "") or "")
        at = str(row.get("assumed_type", "") or "").strip()
        if at == "ET":
            return "L5 ET"
        if at == "Tlx":
            return "L5 IT"  # Tlx3-lineage IT — confirm with user
        if "L2-3" in layer or layer in ("L2/3", "L2-3"):
            return "L2/3 IT"
        if layer == "L5":
            return "L5 (untyped)"
        return "UNKNOWN"

    df["_coarse_type"] = df.apply(map_coarse, axis=1)

    usable = df[~df["exclude_flag"].astype(str).isin(["1", "1.0", "?"])]
    usable = usable[usable["excluded_in_May"].astype(str).isin(["", "nan"]) | usable["excluded_in_May"].isna()]
    print(f"\n--- Usable rows (exclude_flag ok, not excluded_in_May): {len(usable)} / {len(df)} ---")
    print("\nUsable: brain_area_group × coarse_type:")
    print(pd.crosstab(usable["_brain_area_group"], usable["_coarse_type"]).to_string())

    print("\nUsable: brain_area_group × layer:")
    print(pd.crosstab(usable["_brain_area_group"], usable["layer"]).to_string())

    if rf in usable.columns:
        s = pd.to_numeric(usable[rf], errors="coerce")
        print("\nRes freq by brain_area_group × coarse_type (count, median Hz):")
        for (bag, ct), grp in usable.groupby(["_brain_area_group", "_coarse_type"], observed=True):
            vals = pd.to_numeric(grp[rf], errors="coerce").dropna()
            if len(vals):
                print(f"  {bag} × {ct}: n={len(vals)}, median={vals.median():.2f} Hz")

    allen_like = [c for c in df.columns if any(
        x in c.lower() for x in ("supertype", "subclass", "cluster", "allen", "taxonomy")
    )]
    print("\nAllen taxonomy columns:", allen_like or "(none)")

    print("\nDuplicate cell_ids:", int(df["cell_id"].duplicated().sum()))

    rf = "_chirp__Res. freq. (Hz)"
    ri = "_chirp__Res. imp. mag. (MOhm)"
    df["res_freq"] = pd.to_numeric(df[rf], errors="coerce")
    df["res_imp"] = pd.to_numeric(df[ri], errors="coerce")

    print("\n=== Negative or zero res freq ===")
    bad = df[df["res_freq"] <= 0]
    print(f"n={len(bad)}")
    if len(bad):
        cols = ["cell_id", "region", "layer", "assumed_type", "exclude_flag", "res_freq"]
        print(bad[cols].head(20).to_string())

    print("\n=== exclude_flag breakdown ===")
    for flag, grp in df.groupby(df["exclude_flag"].astype(str), dropna=False):
        print(f"  exclude_flag={flag!r}: n={len(grp)}, res_freq median={grp['res_freq'].median():.2f}")

    print("\n=== L5 V2M ET (exclude_flag=0) ===")
    et = df[(df["region"] == "V2M") & (df["layer"] == "L5") & (df["assumed_type"] == "ET")
            & (df["exclude_flag"].astype(str) == "0")]
    print(f"n={len(et)}")
    if len(et):
        print(et[["cell_id", "area_morph", "classic_burster", "res_freq", "res_imp"]].to_string())

    print("\n=== area_morph populated (L5 V2M) ===")
    m = df[(df["region"] == "V2M") & (df["layer"] == "L5") & df["area_morph"].notna()]
    print(m[["cell_id", "area_morph", "area_mismatch", "assumed_type", "res_freq"]].to_string())

    print("\n=== exclude_flag NaN by region×layer ===")
    print(df[df["exclude_flag"].isna()].groupby(["region", "layer"]).size().to_string())


if __name__ == "__main__":
    p = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PATH
    main(p)
