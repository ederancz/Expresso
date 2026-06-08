"""Load and validate receptor_query_config.yaml."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

_REQUIRED_KEYS = ("brain_areas", "cell_type_level", "output", "data")
_VALID_CELL_TYPE_LEVELS = ("class", "subclass", "supertype", "cluster")

# Accepted top-level keys holding the gene panel (family -> [gene symbols]).
# 'gene_panel' is canonical; 'receptors'/'excitability' kept for backward compat.
_GENE_PANEL_KEYS = ("gene_panel", "receptors", "excitability")

# Default root for all notebook outputs (figures, parquet). Outside the git repo.
DEFAULT_OUTPUT_DIR = Path(
    "/Users/rancze/Documents/!Projects/Ach_NE_Marius_Felix/exploration"
)
RUN_MANIFEST_NAME = "run_manifest.json"


def _sanitize_run_slug(label: str) -> str:
    """Filesystem-safe slug for run directory names."""
    slug = "".join(c if c.isalnum() or c in "-_" else "-" for c in label.strip())
    while "--" in slug:
        slug = slug.replace("--", "-")
    return slug.strip("-") or "unknown"


def load_config(path: str | Path = "receptor_query_config.yaml") -> dict[str, Any]:
    """Load YAML config and derive flattened gene lists."""
    config_path = Path(path).expanduser().resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    for key in _REQUIRED_KEYS:
        if key not in cfg:
            raise KeyError(f"Missing required config key: {key}")

    panel_key = next((k for k in _GENE_PANEL_KEYS if k in cfg), None)
    if panel_key is None:
        raise KeyError(
            f"Missing gene panel block; expected one of top-level keys {_GENE_PANEL_KEYS}"
        )
    panel = cfg[panel_key]
    if not isinstance(panel, dict) or not panel:
        raise TypeError(
            f"Gene panel under {panel_key!r} must be a non-empty mapping of "
            "family -> [gene symbols]"
        )

    if cfg["cell_type_level"] not in _VALID_CELL_TYPE_LEVELS:
        raise ValueError(
            f"cell_type_level must be one of {_VALID_CELL_TYPE_LEVELS}, "
            f"got {cfg['cell_type_level']!r}"
        )

    genes: dict[str, str] = {}
    for family, glist in panel.items():
        for g in glist or []:
            genes[g] = family

    cfg["_gene_panel"] = panel
    cfg["_gene_panel_key"] = panel_key
    cfg["_genes_flat"] = genes
    cfg["_all_genes"] = list(genes)
    cfg["_families"] = list(panel.keys())
    cfg["_config_path"] = str(config_path)

    raw_filter = cfg.get("cell_type_name_filter") or []
    if not isinstance(raw_filter, list):
        raise TypeError("cell_type_name_filter must be a list of substrings")
    cfg["cell_type_name_filter"] = [str(s) for s in raw_filter]

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
    panel = config.get("_gene_panel") or config.get("receptors") or {}
    config["_families"] = [
        f for f in config["_families"]
        if any(g in loaded for g in (panel.get(f) or []))
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
    dataset: str,
    exploration_root: Path | str | None = None,
    notebook: str | None = None,
) -> Path:
    """
    Create a timestamped run directory under the exploration root and write
    ``run_manifest.json`` (git repo state + ``receptor_query_config`` snapshot).

    Run folder name: ``{timestamp}_{cell_type_level}_{dataset}``.
    """
    project_root = Path(project_root).resolve()
    root = resolve_output_dir(
        output_dir=exploration_root,
        cfg=cfg if exploration_root is None else None,
    )
    level = cfg["cell_type_level"]
    dataset_slug = _sanitize_run_slug(dataset)
    stamp = (
        datetime.now(timezone.utc).astimezone().strftime("%Y%m%d_%H%M%S")
        + f"_{level}_{dataset_slug}"
    )
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
        "cell_type_level": level,
        "dataset": dataset,
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


DEFAULT_ZHUANG_DATASETS = (
    "Zhuang-ABCA-1",
    "Zhuang-ABCA-2",
    "Zhuang-ABCA-3",
    "Zhuang-ABCA-4",
)


VIZGEN_SAMPLE_TAGS = tuple(f"S{s}R{r}" for s in range(1, 4) for r in range(1, 4))
VIZGEN_FILE_PREFIX = "datasets_mouse_brain_map_BrainReceptorShowcase"


def get_vizgen_data_dir(cfg: dict[str, Any]) -> Path:
    """Return expanded Vizgen flat-file directory."""
    raw = cfg.get("data", {}).get("vizgen_data_dir")
    if not raw:
        raise ValueError(
            "data.vizgen_data_dir is not set; download Vizgen CSVs and set the path "
            "in receptor_query_config.yaml",
        )
    path = Path(raw).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"vizgen_data_dir not found: {path}")
    return path


def vizgen_sample_file_paths(data_dir: Path, sample_tag: str) -> tuple[Path, Path]:
    """
    Return (cell_by_gene, cell_metadata) paths for a Vizgen sample tag (e.g. ``S1R1``).
    """
    import re

    match = re.fullmatch(r"S(\d+)R(\d+)", sample_tag)
    if not match:
        raise ValueError(f"Invalid Vizgen sample tag {sample_tag!r}; expected e.g. 'S1R1'")
    slice_n, rep_n = match.group(1), match.group(2)
    prefix = VIZGEN_FILE_PREFIX
    cbg = (
        data_dir
        / f"{prefix}_Slice{slice_n}_Replicate{rep_n}_cell_by_gene_{sample_tag}.csv"
    )
    meta = (
        data_dir
        / f"{prefix}_Slice{slice_n}_Replicate{rep_n}_cell_metadata_{sample_tag}.csv"
    )
    return cbg, meta


def discover_vizgen_samples(cfg: dict[str, Any]) -> list[str]:
    """List Vizgen sample tags with both cell_by_gene and cell_metadata present."""
    data_dir = get_vizgen_data_dir(cfg)
    found: list[str] = []
    for tag in VIZGEN_SAMPLE_TAGS:
        cbg, meta = vizgen_sample_file_paths(data_dir, tag)
        if cbg.is_file() and meta.is_file():
            found.append(tag)
    return found


def get_vizgen_samples(cfg: dict[str, Any]) -> list[str]:
    """
    Return Vizgen sample tags to process.

    ``data.vizgen_samples``: explicit list (e.g. ``[S1R1, S1R2]``), or ``null`` /
    omitted to auto-discover all downloaded pairs under ``vizgen_data_dir``.
    """
    raw = cfg.get("data", {}).get("vizgen_samples")
    if raw is None:
        samples = discover_vizgen_samples(cfg)
        if not samples:
            raise FileNotFoundError(
                "No Vizgen sample CSV pairs found under "
                f"{get_vizgen_data_dir(cfg)!r}; expected files like "
                f"{VIZGEN_FILE_PREFIX}_Slice1_Replicate1_cell_by_gene_S1R1.csv",
            )
        return samples

    if not isinstance(raw, list) or not raw:
        raise ValueError("data.vizgen_samples must be a non-empty list or null")

    data_dir = get_vizgen_data_dir(cfg)
    samples = [str(s) for s in raw]
    missing: list[str] = []
    for tag in samples:
        cbg, meta = vizgen_sample_file_paths(data_dir, tag)
        if not cbg.is_file() or not meta.is_file():
            missing.append(tag)
    if missing:
        raise FileNotFoundError(
            f"Vizgen samples missing cell_by_gene/cell_metadata: {missing}",
        )
    return samples


def get_zhuang_datasets(cfg: dict[str, Any]) -> list[str]:
    """Return Zhuang replicate dataset IDs from config (default: all four)."""
    raw = cfg.get("data", {}).get("zhuang_datasets")
    if raw is None:
        return list(DEFAULT_ZHUANG_DATASETS)
    if not isinstance(raw, list) or not raw:
        raise ValueError("data.zhuang_datasets must be a non-empty list of dataset IDs")
    return [str(d) for d in raw]


def find_prior_run_parquet(
    cfg: dict[str, Any],
    *,
    parquet_filename: str,
    dataset_slug: str,
    exploration_root: Path | str | None = None,
) -> Path | None:
    """
    Find the most recent prior run parquet under the exploration root.

    Matches run folders named ``{timestamp}_{cell_type_level}_{dataset_slug}``
    (same convention as :func:`start_run`).
    """
    explicit = cfg.get("data", {}).get("allen_merfish_parquet")
    if explicit:
        path = Path(explicit).expanduser()
        if path.is_file():
            return path.resolve()

    root = resolve_output_dir(
        output_dir=exploration_root,
        cfg=cfg if exploration_root is None else None,
    )
    level = cfg["cell_type_level"]
    suffix = f"_{level}_{_sanitize_run_slug(dataset_slug)}"

    candidates: list[Path] = []
    if not root.is_dir():
        return None

    for run_dir in root.iterdir():
        if not run_dir.is_dir() or not run_dir.name.endswith(suffix):
            continue
        parquet = run_dir / parquet_filename
        if parquet.is_file():
            candidates.append(parquet)

    if not candidates:
        return None

    return sorted(candidates, key=lambda p: p.parent.name, reverse=True)[0]
