"""Resolve broad area (VISp/V2M) and areaCCF from messy source strings."""

from __future__ import annotations

import re

from .config import CCF_TO_AREA, VALID_AREA_CCF

_CCF_TOKEN_RE = re.compile(r"(VISpm|VISam|VISp|RSPagl)", re.IGNORECASE)
_BROAD_RE = re.compile(r"\b(V1|V2M|VISp)\b", re.IGNORECASE)


def extract_ccf_tokens(text: str | None) -> list[str]:
    if not text:
        return []
    found: list[str] = []
    for m in _CCF_TOKEN_RE.finditer(str(text)):
        tok = m.group(1)
        # Normalise casing to canonical form
        for valid in VALID_AREA_CCF:
            if tok.lower() == valid.lower():
                tok = valid
                break
        if tok not in found:
            found.append(tok)
    return found


def ccf_to_area(ccf: str) -> str:
    return CCF_TO_AREA.get(ccf, "")


def pick_primary_ccf(tokens: list[str]) -> str:
    """When multiple CCF codes appear, prefer finest VIS subregion over RSPagl."""
    if not tokens:
        return ""
    order = ["VISpm", "VISam", "VISp", "RSPagl"]
    for pref in order:
        if pref in tokens:
            return pref
    return tokens[0]


def parse_morphology_value(raw: str | None) -> tuple[str, str]:
    """Return (area, areaCCF) from a morphology cell value."""
    if raw is None:
        return "", ""
    text = str(raw).strip()
    if not text or text in ("-", "?", "X"):
        return "", ""
    tokens = extract_ccf_tokens(text)
    if not tokens:
        return "", ""
    ccf = pick_primary_ccf(tokens)
    return ccf_to_area(ccf), ccf


def parse_all_cells_area(raw: str | None) -> str:
    if raw is None:
        return ""
    text = str(raw).strip()
    if text.upper() == "V1":
        return "VISp"
    if text.upper() == "V2M":
        return "V2M"
    if text in ("VISp", "V2M"):
        return text
    return ""


def parse_all_cells_note(raw: str | None) -> tuple[str, str]:
    """Return (area, areaCCF) from All cells note column."""
    tokens = extract_ccf_tokens(raw)
    if not tokens:
        return "", ""
    ccf = pick_primary_ccf(tokens)
    return ccf_to_area(ccf), ccf


def resolve_area(
    *,
    sheet_region: str,
    morph_raw: str | None = None,
    all_cells_area: str | None = None,
    all_cells_note: str | None = None,
) -> tuple[str, str, bool]:
    """
    Resolve (area, areaCCF, area_mismatch).

    Priority for areaCCF: morphology row → All cells note → infer from sheet region.
    Priority for broad area: derived from areaCCF → All cells area col → sheet region.
    """
    area_ccf = ""
    area = ""

    morph_area, morph_ccf = parse_morphology_value(morph_raw)
    note_area, note_ccf = parse_all_cells_note(all_cells_note)
    meta_area = parse_all_cells_area(all_cells_area)

    if morph_ccf:
        area_ccf = morph_ccf
        area = morph_area
    elif note_ccf:
        area_ccf = note_ccf
        area = note_area

    if not area:
        if meta_area:
            area = meta_area
        elif sheet_region == "V1":
            area = "VISp"
        elif sheet_region == "V2M":
            area = "V2M"

    if not area_ccf:
        if area == "VISp":
            area_ccf = "VISp"
        elif note_ccf:
            area_ccf = note_ccf

    if not area and area_ccf:
        area = ccf_to_area(area_ccf)

    sources = {x for x in (morph_area, note_area, meta_area) if x}
    mismatch = len(sources) > 1
    return area, area_ccf, mismatch
