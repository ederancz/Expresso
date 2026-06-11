"""Normalise and format electrophysiology parameter values for output."""

from __future__ import annotations

import math
import re
from typing import Any

_NON_NUMERIC_PHRASES = (
    "cannot be determined",
    "rmp cannot",
    "#n/a",
    "#div/0",
    "#ref!",
    "#value!",
)

_EMPTY_MARKERS = frozenset({"", "?", "-", "x"})


def is_non_numeric_marker(val: Any) -> bool:
    if val is None:
        return True
    if isinstance(val, float) and math.isnan(val):
        return True
    if isinstance(val, str):
        s = val.strip()
        if s.lower() in _EMPTY_MARKERS:
            return True
        low = s.lower()
        if any(p in low for p in _NON_NUMERIC_PHRASES):
            return True
        try:
            float(s)
            return False
        except ValueError:
            return True
    return False


def coerce_param_value(val: Any) -> float | None:
    """Return a float or None for parameter columns."""
    if is_non_numeric_marker(val):
        return None
    if isinstance(val, bool):
        return float(int(val))
    if isinstance(val, (int, float)):
        return float(val)
    if isinstance(val, str):
        return float(val.strip())
    return None


def format_sigfigs(val: Any, sigfigs: int = 4) -> Any:
    """Format numeric values to fixed significant figures; None → empty for CSV."""
    if val is None:
        return ""
    if isinstance(val, float) and math.isnan(val):
        return ""
    if isinstance(val, bool):
        return val
    if isinstance(val, int) and not isinstance(val, bool):
        return int(round(float(val), sigfigs)) if abs(val) >= 10 ** (sigfigs - 1) else float(f"{float(val):.{sigfigs}g}")
    if isinstance(val, float):
        return float(f"{val:.{sigfigs}g}")
    return val


def format_param_for_output(val: Any, sigfigs: int = 4) -> Any:
    coerced = coerce_param_value(val)
    if coerced is None:
        return ""
    return format_sigfigs(coerced, sigfigs)
