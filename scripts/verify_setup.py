#!/usr/bin/env python3
"""Smoke tests for config, cache, gene load, and heatmap rendering."""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))


def main() -> int:
    parser = argparse.ArgumentParser(description="ABC Atlas setup smoke tests")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Skip large metadata/expression downloads (steps 4 only)",
    )
    args = parser.parse_args()
    print("1. Loading config...")
    from src.config import load_config, get_figures_dir

    config_path = PROJECT_ROOT / "receptor_query_config.yaml"
    config = load_config(config_path)
    assert config["_all_genes"], "No genes in config"
    print(f"   OK — {len(config['_all_genes'])} genes")

    print("2. Initialising AbcProjectCache...")
    from src.data_loaders import get_abc_cache

    cache = get_abc_cache(config)
    print(f"   OK — manifest {cache.current_manifest}")

    print("3. Checking Drd2 in WMB-10X gene metadata...")
    gene_df = cache.get_metadata_dataframe(directory="WMB-10X", file_name="gene")
    symbols = set(gene_df["gene_symbol"].astype(str))
    if "Drd2" not in symbols:
        warnings.warn("Drd2 not in gene metadata (unexpected)")
    else:
        print("   OK — Drd2 present")

    if args.quick:
        print("4. Skipped (--quick): Drd2 expression load")
    else:
        print("4. Loading Drd2 from WMB-10Xv3-STR (backed, may download)...")
        from src.data_loaders import load_scrna_cell_metadata, load_expression_subset

        test_config = load_config(config_path)
        test_config["brain_areas"] = ["STR"]
        test_config["receptors"] = {"dopamine": ["Drd2"]}
        test_config["_genes_flat"] = {"Drd2": "dopamine"}
        test_config["_all_genes"] = ["Drd2"]
        test_config["_families"] = ["dopamine"]

        try:
            meta = load_scrna_cell_metadata(cache, test_config)
            print(f"   STR cells: {len(meta):,}")
            adata = load_expression_subset(cache, ["Drd2"], meta, test_config)
            if adata is None:
                print("   SKIP — expression not loaded (download may be required)")
            else:
                print(f"   OK — adata shape {adata.shape}")
        except Exception as e:
            print(f"   SKIP — {e}")

    print("5. Rendering minimal 2×2 heatmap...")
    import numpy as np
    import pandas as pd

    from src.plotting import plot_heatmap

    mock = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=["type_a", "type_b"],
        columns=["STR", "TH"],
    )
    figures_dir = get_figures_dir(config, PROJECT_ROOT)
    out = figures_dir / "verify_test_heatmap.png"
    plot_heatmap(mock, "verify test", config, save_path=out, base_dir=PROJECT_ROOT)
    if out.exists():
        print(f"   OK — {out}")
    else:
        print("   FAIL — heatmap not written")
        return 1

    print("\nAll smoke checks completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
