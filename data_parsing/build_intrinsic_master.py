#!/usr/bin/env python3
"""Build Intrinsic_master.xlsx and CSV exports from the Aug2024 analysis workbook."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from intrinsic.build import build_master
from intrinsic.config import DEFAULT_OUTPUT_DIR, DEFAULT_SOURCE
from intrinsic.manifest import (
    print_phase2_pass,
    print_phase3_pass,
    print_phase4_summary,
    print_run_alerts,
)
from intrinsic.verify import IntegrityCheckError
from intrinsic.verify_postbuild import PostBuildVerificationError


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
    try:
        manifest = build_master(args.source, args.output_dir)
    except IntegrityCheckError as exc:
        print(f"\nBUILD ABORTED — {exc.summary()}\n", file=sys.stderr)
        sys.exit(1)
    except PostBuildVerificationError as exc:
        print(f"\nBUILD ABORTED — {exc.summary()}\n", file=sys.stderr)
        sys.exit(1)

    counts = manifest["counts"]

    print(f"control_excitability: {counts['control_neurons']} neurons")
    print(f"pharmacology_effect: {counts['pharmacology_effect_rows']} rows")
    print(f"excluded_in_May dropped: {counts['excluded_in_may_dropped']} cell IDs")
    print(f"parameter columns: {counts['param_columns']}")
    print(f"label merges: {len(manifest['label_merges'])}")
    print(f"cluster fills (assumed_type): {counts['cluster_fills_assumed_type']}")
    print(f"Wrote outputs to {args.output_dir}")
    print(f"run_manifest.json → {args.output_dir / 'run_manifest.json'}")
    print_phase2_pass(manifest)
    print_phase3_pass(manifest)
    print_phase4_summary(manifest)

    print_run_alerts(manifest)


if __name__ == "__main__":
    main()
