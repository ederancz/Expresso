#!/usr/bin/env python3
"""Build Intrinsic_master.xlsx and CSV exports from the Aug2024 analysis workbook."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running as script from repo root or data_parsing/
sys.path.insert(0, str(Path(__file__).resolve().parent))

from intrinsic.build import build_master
from intrinsic.config import DEFAULT_OUTPUT_DIR, DEFAULT_SOURCE


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_SOURCE,
        help="Source Excel workbook",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory for master workbook and CSVs",
    )
    args = parser.parse_args()
    report = build_master(args.source, args.output_dir)
    print(f"control_excitability: {report['control_neurons']} neurons")
    print(f"pharmacology_effect: {report['effect_rows']} rows")
    print(f"duplicate_conflicts: {report['conflict_rows']} rows")
    print(f"parameter columns: {report['param_columns']}")
    print(f"label merges: {len(report['label_merges'])}")
    print(f"cluster fills (assumed_type): {report['cluster_fills_assumed_type']}")
    print(f"Wrote outputs to {args.output_dir}")


if __name__ == "__main__":
    main()
