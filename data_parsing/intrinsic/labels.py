"""Parameter label union and fuzzy merge."""

from __future__ import annotations

import re
from difflib import SequenceMatcher

from .config import EXPLICIT_LABEL_MERGES, LABEL_MERGE_RATIO, SKIP_PARAM_LABELS

_DIGIT_RE = re.compile(r"\d+")


def normalize_section(header: str) -> str:
    h = header.strip()
    if "_IV" in h:
        return "_IV"
    if h.startswith("_short"):
        return "_short_depol"
    if h.startswith("_EPSP"):
        return "_EPSP"
    if h.startswith("_crit"):
        return "_crit_freq"
    if h.startswith("_sag"):
        return "_sag"
    if h.startswith("_chirp"):
        return "_chirp"
    return h.split()[0]


def iv_header_label(header: str) -> str | None:
    if "Vm_from_fit" in header:
        return "Vm_from_fit (mV)"
    return None


def col_name(section: str, label: str) -> str:
    safe = label.replace(" ", "_")
    return f"{section}__{safe}"


def _should_fuzzy_merge(a: str, b: str, ratio: float) -> bool:
    if ratio < LABEL_MERGE_RATIO:
        return False
    # Distinct numbered variants (e.g. thresh vs thresh2) must not merge.
    if _DIGIT_RE.findall(a) != _DIGIT_RE.findall(b):
        return False
    return True


def build_label_merge_map(raw_pairs: set[tuple[str, str]]) -> tuple[dict[tuple[str, str], str], list[str]]:
    """Return (merge_map, log_lines). merge_map maps (section, raw) → canonical."""
    merge_map: dict[tuple[str, str], str] = {}
    log: list[str] = []

    for section, raw in sorted(raw_pairs):
        if (section, raw) in EXPLICIT_LABEL_MERGES:
            canonical = EXPLICIT_LABEL_MERGES[(section, raw)]
            if canonical != raw:
                merge_map[(section, raw)] = canonical
                log.append(f"explicit: ({section}, {raw!r}) → {canonical!r}")

    by_section: dict[str, list[str]] = {}
    for section, raw in raw_pairs:
        canonical = merge_map.get((section, raw), raw)
        by_section.setdefault(section, [])
        if canonical not in by_section[section]:
            by_section[section].append(canonical)

    for section, labels in by_section.items():
        labels = sorted(set(labels))
        for i, a in enumerate(labels):
            for b in labels[i + 1 :]:
                ratio = SequenceMatcher(None, a, b).ratio()
                if _should_fuzzy_merge(a, b, ratio) and a != b:
                    # Prefer longer / cleaner label as canonical
                    canonical, variant = (a, b) if len(a) >= len(b) else (b, a)
                    if (section, variant) not in merge_map:
                        merge_map[(section, variant)] = canonical
                        log.append(f"fuzzy ({ratio:.3f}): ({section}, {variant!r}) → {canonical!r}")

    return merge_map, log


def canonical_label(section: str, raw: str, merge_map: dict[tuple[str, str], str]) -> str | None:
    if raw in SKIP_PARAM_LABELS:
        return None
    return merge_map.get((section, raw), raw)
