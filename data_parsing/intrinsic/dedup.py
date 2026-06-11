"""Deduplication at 4 significant figures."""

from __future__ import annotations

from typing import Any

from .config import DEDUP_PRIORITY


def _priority_index(source_sheet: str) -> int:
    try:
        return DEDUP_PRIORITY.index(source_sheet)
    except ValueError:
        return len(DEDUP_PRIORITY)


def values_equal(a: Any, b: Any, sigfigs: int = 4) -> bool:
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    if type(a) != type(b):
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            pass
        else:
            return str(a) == str(b)
    if isinstance(a, bool):
        return a == b
    if isinstance(a, str):
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a == b:
            return True
        if a == 0 or b == 0:
            return abs(a - b) < 10 ** (-sigfigs)
        scale = max(abs(a), abs(b))
        if scale == 0:
            return True
        rel = abs(a - b) / scale
        return rel < 10 ** (-sigfigs)
    return a == b


def params_agree(
    params_a: dict[str, Any],
    params_b: dict[str, Any],
    sigfigs: int = 4,
) -> bool:
    keys = set(params_a) & set(params_b)
    for k in keys:
        if not values_equal(params_a[k], params_b[k], sigfigs):
            return False
    return True


def conflicting_params(
    params_a: dict[str, Any],
    params_b: dict[str, Any],
    source_a: str,
    source_b: str,
    sigfigs: int = 4,
) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for k in set(params_a) & set(params_b):
        if not values_equal(params_a[k], params_b[k], sigfigs):
            out[k] = {source_a: params_a[k], source_b: params_b[k]}
    return out


def pick_canonical(instances: list[dict[str, Any]]) -> dict[str, Any]:
    """Pick canonical instance by source priority; merge non-overlapping keys."""
    ordered = sorted(instances, key=lambda x: _priority_index(x["source_sheet"]))
    canonical = dict(ordered[0])
    params = dict(canonical.get("params", {}))
    for inst in ordered[1:]:
        for k, v in inst.get("params", {}).items():
            if k not in params:
                params[k] = v
    canonical["params"] = params
    canonical["sheet_meta"] = merge_sheet_metas(instances)
    for inst in ordered:
        if inst.get("region"):
            canonical["region"] = inst["region"]
            canonical["layer"] = inst["layer"]
            break
    return canonical


def merge_sheet_metas(instances: list[dict[str, Any]]) -> dict[str, dict]:
    classic: dict[str, Any] = {}
    morph: dict[str, Any] = {}
    for inst in instances:
        sm = inst.get("sheet_meta", {})
        for cid, v in sm.get("classic_burster", {}).items():
            if v is not None and cid not in classic:
                classic[cid] = v
        for cid, v in sm.get("area_morph_raw", {}).items():
            if v is not None and str(v).strip() not in ("", "-", "?", "X"):
                morph[cid] = v
    return {"classic_burster": classic, "area_morph_raw": morph}


def deduplicate_control_instances(
    instances: list[dict[str, Any]],
) -> tuple[dict[str, Any] | None, list[dict[str, Any]], bool, list[dict[str, Any]]]:
    """
    Returns (canonical, conflict_instances, has_conflict, conflict_detail).
    conflict_detail entries: {cell_id, sources, conflicting_params}.
    """
    if not instances:
        return None, [], False, []
    if len(instances) == 1:
        return instances[0], [], False, []

    detail: list[dict[str, Any]] = []
    has_conflict = False
    for i, a in enumerate(instances):
        for b in instances[i + 1 :]:
            disagree = conflicting_params(
                a.get("params", {}),
                b.get("params", {}),
                a["source_sheet"],
                b["source_sheet"],
            )
            if disagree:
                has_conflict = True
                detail.append(
                    {
                        "sources": [a["source_sheet"], b["source_sheet"]],
                        "conflicting_params": disagree,
                    }
                )

    canonical = pick_canonical(instances)
    if has_conflict:
        return canonical, list(instances), True, detail
    return canonical, [], False, []


def summarize_conflicts(
    all_details: list[dict[str, Any]],
) -> tuple[int, int]:
    """Return (n_cells, n_instances) with measured-value conflicts."""
    cells = {d["cell_id"] for d in all_details}
    instances = sum(len(d.get("instances", [])) for d in all_details)
    return len(cells), instances
