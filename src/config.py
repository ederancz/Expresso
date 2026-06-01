"""Load and validate receptor_query_config.yaml."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

_REQUIRED_KEYS = ("receptors", "brain_areas", "cell_type_level", "output", "data")
_VALID_CELL_TYPE_LEVELS = ("class", "subclass", "supertype", "cluster")


def load_config(path: str | Path = "receptor_query_config.yaml") -> dict[str, Any]:
    """Load YAML config and derive flattened gene lists."""
    config_path = Path(path).expanduser().resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    for key in _REQUIRED_KEYS:
        if key not in cfg:
            raise KeyError(f"Missing required config key: {key}")

    if cfg["cell_type_level"] not in _VALID_CELL_TYPE_LEVELS:
        raise ValueError(
            f"cell_type_level must be one of {_VALID_CELL_TYPE_LEVELS}, "
            f"got {cfg['cell_type_level']!r}"
        )

    genes: dict[str, str] = {}
    for family, glist in cfg["receptors"].items():
        for g in glist:
            genes[g] = family

    cfg["_genes_flat"] = genes
    cfg["_all_genes"] = list(genes)
    cfg["_families"] = list(cfg["receptors"].keys())
    cfg["_config_path"] = str(config_path)

    return cfg


def get_figures_dir(cfg: dict[str, Any], base_dir: Path | None = None) -> Path:
    """Return figures output directory, creating it if needed."""
    rel = cfg["output"]["figures_dir"]
    root = base_dir or Path(cfg["_config_path"]).parent
    out = (root / rel).resolve()
    out.mkdir(parents=True, exist_ok=True)
    return out


def get_parquet_path(cfg: dict[str, Any], base_dir: Path | None = None) -> Path:
    """Return path for aggregated scRNA parquet cache."""
    root = base_dir or Path(cfg["_config_path"]).parent
    data_dir = root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir / "aggregated_scrna.parquet"


def get_cache_dir(cfg: dict[str, Any]) -> Path:
    """Return expanded ABC Atlas cache directory."""
    return Path(cfg["data"]["cache_dir"]).expanduser().resolve()


def get_expression_suffix(cfg: dict[str, Any]) -> str:
    """Return 'log2' or 'raw' for expression matrix file paths."""
    unit = cfg["data"].get("expression_unit", "log2")
    if unit not in ("log2", "raw"):
        raise ValueError(f"expression_unit must be 'log2' or 'raw', got {unit!r}")
    return unit
