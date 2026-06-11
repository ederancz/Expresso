"""Run manifest for intrinsic master builds."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import REPO_ROOT

RUN_MANIFEST_NAME = "run_manifest.json"


def collect_git_info() -> dict[str, Any]:
    root = REPO_ROOT.resolve()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    try:
        from src.config import collect_git_info as _collect

        return _collect(root)
    except Exception:
        return {"root": str(root), "commit": None, "commit_short": None, "branch": None}


def write_run_manifest(
    output_dir: Path,
    *,
    source_path: Path,
    report: dict[str, Any],
) -> dict[str, Any]:
    """Write ``run_manifest.json`` and return the manifest dict."""
    manifest: dict[str, Any] = {
        "built_at": datetime.now(timezone.utc).astimezone().isoformat(),
        "git": collect_git_info(),
        "source_workbook": str(source_path.resolve()),
        "outputs": {
            "intrinsic_master_xlsx": str((output_dir / "Intrinsic_master.xlsx").resolve()),
            "control_excitability_csv": str((output_dir / "control_excitability.csv").resolve()),
            "pharmacology_effect_csv": str((output_dir / "pharmacology_effect.csv").resolve()),
            "duplicate_conflicts_csv": str((output_dir / "duplicate_conflicts.csv").resolve()),
        },
        "counts": {
            "control_neurons": report.get("control_neurons", 0),
            "pharmacology_effect_rows": report.get("effect_rows", 0),
            "duplicate_conflict_instances": report.get("conflict_rows", 0),
            "duplicate_conflict_cells": report.get("conflict_cells", 0),
            "excluded_in_may_dropped": report.get("excluded_in_may_dropped", 0),
            "param_columns": report.get("param_columns", 0),
            "cluster_fills_assumed_type": report.get("cluster_fills_assumed_type", 0),
        },
        "cluster_fill_column": report.get("cluster_fill_column", "assumed_type"),
        "label_merges": report.get("label_merges", []),
        "duplicate_header_warnings": report.get("duplicate_header_warnings", []),
        "duplicate_conflicts": report.get("duplicate_conflicts_detail", []),
        "region_area_conflicts": report.get("region_area_conflicts", []),
        "notes": {
            "duplicate_header_warnings": (
                "CRITICAL: the same cell ID appears as two column headers on one sheet "
                "(likely a student typo). Parsed as separate IDs with #1/#2 suffixes. "
                "Review before analysis."
            ),
            "dup_conflict": (
                "dup_conflict=TRUE when the same cell_id appears on multiple source sheets "
                "with overlapping parameter values that disagree beyond 4 significant figures. "
                "Review duplicate_conflicts.csv and duplicate_conflicts in this manifest."
            ),
            "exclude_flag": (
                "Informational only — copied verbatim from the source workbook. "
                "Not used to exclude cells in downstream Expresso analysis."
            ),
            "region": (
                "Broad atlas region (VISp or V2M). V1 from source sheets is normalised to VISp."
            ),
        },
    }
    path = output_dir / RUN_MANIFEST_NAME
    path.write_text(json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8")
    return manifest


def print_run_alerts(manifest: dict[str, Any]) -> None:
    """Print build alerts to terminal — duplicate headers loudest, then measurement conflicts."""
    banner = "=" * 78
    header_warnings = manifest.get("duplicate_header_warnings") or []
    counts = manifest.get("counts", {})
    n_meas_cells = counts.get("duplicate_conflict_cells", 0)
    n_meas_inst = counts.get("duplicate_conflict_instances", 0)
    region_conflicts = manifest.get("region_area_conflicts") or []

    if header_warnings:
        print(f"\n{banner}", flush=True)
        print(
            "*** CRITICAL: DUPLICATE CELL ID HEADERS IN SOURCE WORKBOOK ***",
            flush=True,
        )
        print(
            "    Same cell ID in two columns — almost certainly a student typo.",
            flush=True,
        )
        print(
            "    Parsed as separate neurons (#1 / #2 suffixes); verify before trusting data.",
            flush=True,
        )
        print(f"{banner}", flush=True)
        for w in header_warnings:
            print(f"  • {w}", flush=True)
        print(
            f"\n  → {len(header_warnings)} warning(s): run_manifest.json → duplicate_header_warnings\n",
            flush=True,
        )

    if n_meas_cells:
        print(f"{banner}", flush=True)
        print(
            f"*** DUPLICATE MEASUREMENT CONFLICTS (RED FLAG): "
            f"{n_meas_cells} cell(s), {n_meas_inst} conflicting instance row(s) ***",
            flush=True,
        )
        print(
            "    Same cell_id on multiple sheets with disagreeing parameter values (>4 sig figs).",
            flush=True,
        )
        print(
            "  → duplicate_conflicts.csv  |  run_manifest.json → duplicate_conflicts",
            flush=True,
        )
        print(f"{banner}\n", flush=True)
    elif not header_warnings:
        print("\nNo duplicate measurement conflicts detected.", flush=True)

    if region_conflicts:
        print(f"{banner}", flush=True)
        print(
            f"*** REGION CONFLICTS: {len(region_conflicts)} cell(s) "
            f"where sheet region ≠ metadata/CCF-derived region ***",
            flush=True,
        )
        print("  → run_manifest.json → region_area_conflicts", flush=True)
        print(f"{banner}\n", flush=True)
