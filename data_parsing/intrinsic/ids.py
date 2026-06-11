"""Cell ID normalisation and TASK disambiguation."""

from __future__ import annotations

import re
from dataclasses import dataclass, field

_ID_RE = re.compile(r"^(nm)?(\d{4}_\d{2}_\d{2}_c\d+)", re.IGNORECASE)


@dataclass
class HeaderInfo:
    col: int
    raw: str
    normalized: str
    task_note: str = ""


def normalize_id(raw: str) -> str:
    s = str(raw).strip()
    m = _ID_RE.match(s)
    if not m:
        return s
    body = m.group(2)
    return f"nm{body}"


def parse_headers_row(values: list[tuple[int, object]], *, disambiguate_dupes: bool = False) -> list[HeaderInfo]:
    """Parse (col_index, header_value) pairs into HeaderInfo list."""
    parsed: list[tuple[int, str, str]] = []
    for col, raw in values:
        if raw is None:
            continue
        text = str(raw).strip()
        if not _ID_RE.match(text) and not text.lower().startswith("nm"):
            continue
        parsed.append((col, text, normalize_id(text)))

    if disambiguate_dupes:
        counts: dict[str, int] = {}
        for _, _, norm in parsed:
            counts[norm] = counts.get(norm, 0) + 1
        seen: dict[str, int] = {}
        out: list[HeaderInfo] = []
        for col, text, norm in parsed:
            note = ""
            if counts[norm] > 1:
                seen[norm] = seen.get(norm, 0) + 1
                suffix = seen[norm]
                note = f"duplicate header disambiguated as #{suffix}"
                norm = f"{norm}#{suffix}"
            out.append(HeaderInfo(col=col, raw=text, normalized=norm, task_note=note))
        return out

    return [HeaderInfo(col=col, raw=text, normalized=norm) for col, text, norm in parsed]


def duplicate_header_warnings(headers: list[HeaderInfo]) -> list[str]:
    """Prominent warnings for disambiguated duplicate headers (likely student typos)."""
    warnings: list[str] = []
    for h in headers:
        if "#" in h.normalized:
            warnings.append(
                f"DUPLICATE CELL ID HEADER: raw={h.raw!r} → {h.normalized!r} "
                f"(same ID in two columns — likely a student typo; review source sheet)"
            )
    return warnings
