"""Load and validate receptor_query_config.yaml."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

_REQUIRED_KEYS = ("receptors", "brain_areas", "cell_type_level", "output", "data")
_VALID_CELL_TYPE_LEVELS = ("class", "subclass", "supertype", "cluster")

# Default root for all notebook outputs (figures, parquet). Outside the git repo.
DEFAULT_OUTPUT_DIR = Path(
    "/Users/rancze/Documents/!Projects/Ach_NE_Marius_Felix/exploration"
)
RUN_MANIFEST_NAME = "run_manifest.json"


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


def restrict_config_to_genes(config: dict[str, Any], gene_symbols: list[str]) -> list[str]:
    """
    Restrict derived config gene lists to ``gene_symbols``.

    Returns sorted list of previously requested symbols that were removed.
    """
    loaded = set(gene_symbols)
    requested = list(config["_all_genes"])
    missing = sorted(set(requested) - loaded)

    config["_all_genes"] = [g for g in requested if g in loaded]
    config["_genes_flat"] = {
        g: f for g, f in config["_genes_flat"].items() if g in loaded
    }
    config["_families"] = [
        f for f in config["_families"]
        if any(g in loaded for g in config["receptors"].get(f, []))
    ]
    return missing


def _git_command(repo_root: Path, *args: str) -> str | None:
    try:
        return (
            subprocess.check_output(
                ["git", "-C", str(repo_root), *args],
                text=True,
                stderr=subprocess.DEVNULL,
            )
            .strip()
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def collect_git_info(project_root: Path | str) -> dict[str, Any]:
    """Collect git metadata for ``run_manifest.json``."""
    root = Path(project_root).resolve()
    commit = _git_command(root, "rev-parse", "HEAD")
    dirty_output = _git_command(root, "status", "--porcelain")
    return {
        "root": str(root),
        "remote_url": _git_command(root, "remote", "get-url", "origin"),
        "branch": _git_command(root, "rev-parse", "--abbrev-ref", "HEAD"),
        "commit": commit,
        "commit_short": _git_command(root, "rev-parse", "--short", "HEAD"),
        "dirty": bool(dirty_output),
    }


def start_run(
    project_root: Path | str,
    cfg: dict[str, Any],
    *,
    exploration_root: Path | str | None = None,
    notebook: str | None = None,
) -> Path:
    """
    Create a timestamped run directory under the exploration root and write
    ``run_manifest.json`` (git repo state + ``receptor_query_config`` snapshot).
    """
    project_root = Path(project_root).resolve()
    root = resolve_output_dir(
        output_dir=exploration_root,
        cfg=cfg if exploration_root is None else None,
    )
    stamp = datetime.now(timezone.utc).astimezone().strftime("%Y%m%d_%H%M%S")
    run_dir = root / stamp
    suffix = 0
    while run_dir.exists():
        suffix += 1
        run_dir = root / f"{stamp}_{suffix}"

    config_path = Path(cfg["_config_path"])
    with open(config_path) as f:
        config_contents = yaml.safe_load(f)

    manifest: dict[str, Any] = {
        "run_id": run_dir.name,
        "started_at": datetime.now(timezone.utc).astimezone().isoformat(),
        "notebook": notebook,
        "exploration_root": str(root),
        "run_dir": str(run_dir.resolve()),
        "git": collect_git_info(project_root),
        "receptor_query_config": {
            "path": str(config_path),
            "contents": config_contents,
        },
    }
    run_dir.mkdir(parents=True)
    with open(run_dir / RUN_MANIFEST_NAME, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return run_dir.resolve()


def resolve_output_dir(
    output_dir: Path | str | None = None,
    cfg: dict[str, Any] | None = None,
    *,
    base_dir: Path | None = None,
) -> Path:
    """Return root directory for all notebook outputs, creating it if needed.

    Resolution order: explicit ``output_dir`` → ``output.output_dir`` in YAML
    (absolute paths as-is; relative paths under ``base_dir`` or config parent)
    → :data:`DEFAULT_OUTPUT_DIR`.
    """
    if output_dir is not None:
        out = Path(output_dir).expanduser().resolve()
    elif cfg is not None:
        raw = cfg.get("output", {}).get("output_dir") or cfg.get("output", {}).get(
            "figures_dir"
        )
        if raw:
            path = Path(raw).expanduser()
            if path.is_absolute():
                out = path.resolve()
            else:
                root = base_dir or Path(cfg["_config_path"]).parent
                out = (root / path).resolve()
        else:
            out = DEFAULT_OUTPUT_DIR.resolve()
    else:
        out = DEFAULT_OUTPUT_DIR.resolve()
    out.mkdir(parents=True, exist_ok=True)
    return out


def get_figures_dir(
    cfg: dict[str, Any],
    base_dir: Path | None = None,
    output_dir: Path | str | None = None,
) -> Path:
    """Return figures output directory (same root as other notebook outputs)."""
    return resolve_output_dir(
        output_dir=output_dir,
        cfg=cfg if output_dir is None else None,
        base_dir=base_dir,
    )


def get_parquet_path(
    cfg: dict[str, Any],
    output_dir: Path | str | None = None,
    *,
    base_dir: Path | None = None,
) -> Path:
    """Return path for aggregated scRNA parquet cache (under output dir, not repo)."""
    root = resolve_output_dir(
        output_dir=output_dir,
        cfg=cfg if output_dir is None else None,
        base_dir=base_dir,
    )
    return root / "aggregated_scrna.parquet"


def get_cache_dir(cfg: dict[str, Any]) -> Path:
    """Return expanded ABC Atlas cache directory."""
    return Path(cfg["data"]["cache_dir"]).expanduser().resolve()


def get_expression_suffix(cfg: dict[str, Any]) -> str:
    """Return 'log2' or 'raw' for expression matrix file paths."""
    unit = cfg["data"].get("expression_unit", "log2")
    if unit not in ("log2", "raw"):
        raise ValueError(f"expression_unit must be 'log2' or 'raw', got {unit!r}")
    return unit
