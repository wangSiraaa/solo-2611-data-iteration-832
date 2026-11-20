"""Exact closed rational arithmetic shared by the ledger and the quantitation vault.

Everything downstream (concentrations, masses, LP vertices) is a :class:`Fraction`.
Inputs are parsed from decimal or ``a/b`` strings, never from binary floats.
"""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
from typing import Iterable

Q = Fraction

ZERO = Fraction(0)
ONE = Fraction(1)


def q(value, field: str = "value") -> Fraction:
    """Parse a positive/zero rational from int, str or Decimal.

    Floats are rejected on purpose: ``0.1`` is not the rational 1/10 and would
    poison the exact endpoint guarantee.
    """
    if isinstance(value, bool):
        raise ValueError(f"{field} must be a rational number, not a bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, Decimal):
        return Fraction(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError(f"{field} is empty")
        try:
            return Fraction(text)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{field}={value!r} is not a rational (use '1/3' or '0.25')") from exc
    raise ValueError(f"{field} must be an int, str, Decimal or Fraction, got {type(value).__name__}")


def qstr(value: Fraction) -> str:
    """Canonical exact spelling: integers as ``3``, otherwise ``a/b``."""
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def qjson(value: Fraction) -> str:
    return qstr(value)


def nonneg(value: Fraction, field: str) -> Fraction:
    if value < 0:
        raise ValueError(f"{field} must be non-negative, got {qstr(value)}")
    return value


def positive(value: Fraction, field: str) -> Fraction:
    if value <= 0:
        raise ValueError(f"{field} must be strictly positive, got {qstr(value)}")
    return value


def fmt_vec(vec: Iterable[Fraction]) -> str:
    return "(" + ", ".join(qstr(v) for v in vec) + ")"
