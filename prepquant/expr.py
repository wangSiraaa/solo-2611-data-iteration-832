"""Dense exact linear expressions in the sample unknowns.

A container's volume is a constant (volumes never depend on concentration).
A container's solute mass is ``const + sum(coeff_i * x_i)`` where ``x_i`` is the
unknown concentration of sample ``i`` (mg/mL). Unknowns are identified by their
name; the dense coefficient order is fixed by the ledger snapshot.

Mass instead of concentration is the natural currency because every operation in
an ideal solution (transfer, mix, evaporate, top-up) is linear in mass and the
concentration is recovered as ``mass / volume`` exactly.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Mapping, Sequence

from .rational import ZERO, qstr


class LinExpr:
    __slots__ = ("const", "coeffs", "unknowns")

    def __init__(
        self,
        const: Fraction = ZERO,
        coeffs: Sequence[Fraction] = (),
        unknowns: tuple[str, ...] = (),
    ) -> None:
        self.const = Fraction(const)
        self.coeffs = tuple(Fraction(c) for c in coeffs)
        self.unknowns = tuple(unknowns)
        if self.coeffs and len(self.coeffs) != len(self.unknowns):
            raise ValueError("coefficient vector does not match unknown names")

    # -- construction -----------------------------------------------------
    @classmethod
    def constant(cls, value: Fraction, unknowns: Iterable[str]) -> "LinExpr":
        names = tuple(unknowns)
        return cls(Fraction(value), (ZERO,) * len(names), names)

    @classmethod
    def variable(cls, name: str, unknowns: Iterable[str], scale: Fraction = Fraction(1)) -> "LinExpr":
        names = tuple(unknowns)
        idx = names.index(name)
        coeffs = [ZERO] * len(names)
        coeffs[idx] = Fraction(scale)
        return cls(ZERO, coeffs, names)

    @classmethod
    def from_rows(
        cls,
        const: Fraction,
        coeffs: Mapping[str, Fraction],
        unknowns: Iterable[str],
    ) -> "LinExpr":
        names = tuple(unknowns)
        vec = [Fraction(coeffs.get(name, ZERO)) for name in names]
        return cls(const, vec, names)

    def _align(self, other: "LinExpr") -> tuple["LinExpr", "LinExpr"]:
        if not self.unknowns:
            return LinExpr(self.const, (ZERO,) * len(other.unknowns), other.unknowns), other
        if not other.unknowns:
            return self, LinExpr(other.const, (ZERO,) * len(self.unknowns), self.unknowns)
        if self.unknowns != other.unknowns:
            raise ValueError("expressions refer to different unknown sets")
        return self, other

    # -- arithmetic -------------------------------------------------------
    def __add__(self, other: "LinExpr") -> "LinExpr":
        a, b = self._align(other)
        return LinExpr(
            a.const + b.const,
            [x + y for x, y in zip(a.coeffs, b.coeffs)],
            a.unknowns,
        )

    def __sub__(self, other: "LinExpr") -> "LinExpr":
        a, b = self._align(other)
        return LinExpr(
            a.const - b.const,
            [x - y for x, y in zip(a.coeffs, b.coeffs)],
            a.unknowns,
        )

    def scale(self, factor: Fraction) -> "LinExpr":
        factor = Fraction(factor)
        return LinExpr(self.const * factor, [c * factor for c in self.coeffs], self.unknowns)

    def shift(self, delta: Fraction) -> "LinExpr":
        return LinExpr(self.const + Fraction(delta), self.coeffs, self.unknowns)

    # -- queries ----------------------------------------------------------
    def evaluate(self, values: Mapping[str, Fraction]) -> Fraction:
        total = self.const
        for name, coeff in zip(self.unknowns, self.coeffs):
            if coeff:
                total += coeff * values[name]
        return total

    def coeff_of(self, name: str) -> Fraction:
        if name not in self.unknowns:
            return ZERO
        return self.coeffs[self.unknowns.index(name)]

    def rebase(self, unknowns: Iterable[str]) -> "LinExpr":
        """Re-express over a (possibly wider) ordered unknown set."""
        names = tuple(unknowns)
        values = {n: c for n, c in zip(self.unknowns, self.coeffs)}
        return LinExpr(self.const, [values.get(n, ZERO) for n in names], names)

    def is_constant(self) -> bool:
        return all(c == 0 for c in self.coeffs)

    def has_unknown(self) -> bool:
        return not self.is_constant()

    def sparse(self) -> dict[str, Fraction]:
        return {n: c for n, c in zip(self.unknowns, self.coeffs) if c != 0}

    def formatted(self, variable: str | None = None) -> str:
        """e.g. ``12 + 5/2*s1 - 1/3*s2`` (leading zero terms dropped).

        ``variable`` may be given to render positional names (``x1``, ``x2``)
        instead of the declared sample names.
        """
        parts: list[str] = []
        if self.const != 0 or self.is_constant():
            parts.append(qstr(self.const))
        for i, coeff in enumerate(self.coeffs):
            if coeff == 0:
                continue
            name = f"{variable}{i + 1}" if variable else self.unknowns[i]
            term = _coef_term(coeff, name)
            if not parts:
                parts.append(term)
            elif term.startswith("-"):
                parts.append("- " + term[1:])
            else:
                parts.append("+ " + term)
        return " ".join(parts) if parts else "0"

    def __repr__(self) -> str:  # pragma: no cover - debug only
        return f"LinExpr({self.formatted()})"


def _coef_term(coeff: Fraction, var: str) -> str:
    if coeff == 1:
        return var
    if coeff == -1:
        return "-" + var
    return f"{qstr(coeff)}*{var}"
