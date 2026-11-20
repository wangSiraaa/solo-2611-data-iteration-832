"""Exact two-phase simplex over :class:`fractions.Fraction`.

No floating point is used anywhere, so every vertex and every feasibility verdict
is exact. The polytope dimensions are tiny (at most six samples, twenty
measurement rows plus bound rows), so a dense tableau with ``Bland``'s rule is
both fast and provably free of cycling.

Standard form: ``M z = b``, ``z >= 0``. Structural variables are the sample
concentrations (all non-negative by the problem statement). Bounds are added as
explicit inequality rows; equalities become equality rows with an artificial
variable that phase I drives to zero.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

from .rational import ZERO, qstr


class Infeasible(Exception):
    pass


class Unbounded(Exception):
    pass


class SolverError(Exception):
    pass


INF = float("inf")


def _new_col(rows: list[dict[int, Fraction]]) -> int:
    """Return a fresh column index (all rows are dict-sparse until finalised)."""
    width = 0
    for row in rows:
        if row:
            width = max(width, max(row) + 1)
    return width


def build_standard_form(
    rows: Sequence[tuple[Sequence[Fraction], str, Fraction]],
    n: int,
    bounds: Sequence[tuple[Fraction | None, Fraction | None]],
) -> tuple[list[list[Fraction]], list[Fraction], list[int], set[int], int]:
    """Translate constraints into ``M z = b`` with **all** slacks/surplus present.

    Columns ``0..n-1`` are the sample concentrations; every comparison row then
    gets its own slack/surplus column (allocated above ``n``). Equality rows have
    no surplus. This gives one sided ``b >= 0`` whenever the original right-hand
    side is positive, so the all-slack basis is feasible without artificials.
    Returns ``(M, b, basis, artificial_cols, width)`` (``artificial`` is empty
    unless an interval bound forces a negative right-hand side).
    """
    sparse_rows: list[dict[int, Fraction]] = []
    rhs: list[Fraction] = []
    basis: list[int] = []
    artificial: set[int] = set()
    next_col = n

    def new_col() -> int:
        nonlocal next_col
        col = next_col
        next_col += 1
        return col

    def add_row(vec: dict[int, Fraction], b: Fraction, natural: int | None) -> None:
        if b < 0:
            vec = {j: -v for j, v in vec.items()}
            b = -b
            natural = None if natural is not None and natural not in vec else natural
        if natural is not None and vec.get(natural) == 1:
            basis.append(natural)
        else:
            art = new_col()
            vec = dict(vec)
            vec[art] = Fraction(1)
            artificial.add(art)
            basis.append(art)
        sparse_rows.append(vec)
        rhs.append(b)

    def add_constraint(coeffs: Sequence[Fraction], op: str, r: Fraction) -> None:
        vec = {j: Fraction(c) for j, c in enumerate(coeffs) if c != 0}
        if op == "<=":
            slack = new_col()
            vec[slack] = Fraction(1)
            add_row(vec, Fraction(r), slack)
        elif op == ">=":
            surplus = new_col()
            vec[surplus] = Fraction(-1)
            add_row(vec, Fraction(r), None)
        elif op == "=":
            add_row(vec, Fraction(r), None)
        else:  # pragma: no cover
            raise ValueError(f"unknown comparison {op!r}")

    for coeffs, op, r in rows:
        add_constraint(coeffs, op, Fraction(r))

    # Upper bounds get a slack column; lower bounds above the structural x >= 0
    # become ``x >= lo`` with a surplus column (an artificial appears only for
    # the lower-bound rows, and phase I drives them out).
    for j, (lo, hi) in enumerate(bounds):
        if lo is not None and lo > 0:
            add_constraint([Fraction(1) if k == j else ZERO for k in range(n)], ">=", Fraction(lo))
        if hi is not None:
            vec = {j: Fraction(1)}
            slack = new_col()
            vec[slack] = Fraction(1)
            add_row(vec, Fraction(hi), slack)

    width = next_col
    M = [[row.get(j, ZERO) for j in range(width)] for row in sparse_rows]
    return M, rhs, basis, artificial, width


def _matvec(Binv: list[list[Fraction]], vec: list[Fraction]) -> list[Fraction]:
    return [sum((a * b for a, b in zip(row, vec)), ZERO) for row in Binv]


def _row_canonical(Binv: list[list[Fraction]], M: list[list[Fraction]], r: int) -> list[Fraction]:
    width = len(M[0]) if M else 0
    out = []
    for j in range(width):
        out.append(sum((Binv[r][i] * M[i][j] for i in range(len(M))), ZERO))
    return out


def _reduced_costs(
    M: list[list[Fraction]],
    Binv: list[list[Fraction]],
    c: list[Fraction],
    basis: list[int],
) -> list[Fraction]:
    m = len(M)
    width = len(c)
    y = [ZERO] * m
    for k in range(m):
        ck = c[basis[k]]
        if ck:
            brow = Binv[k]
            for i in range(m):
                if brow[i]:
                    y[i] += ck * brow[i]
    rc = []
    for j in range(width):
        val = c[j]
        for i in range(m):
            if y[i] and M[i][j]:
                val -= y[i] * M[i][j]
        rc.append(val)
    return rc


def _pivot(
    Binv: list[list[Fraction]],
    basis: list[int],
    row: int,
    col: int,
    direction: list[Fraction],
) -> None:
    """Replace ``basis[row]`` by ``col`` and update ``Binv``.

    ``direction`` is ``Binv @ A_col``; its ``row`` entry is the pivot. The row
    operations eliminate the *direction* components, not raw ``Binv`` entries.
    """
    pivot = direction[row]
    if pivot == 0:
        raise SolverError("zero pivot")
    m = len(basis)
    for k in range(m):
        if k == row:
            continue
        factor = direction[k] / pivot
        if factor:
            src = Binv[row]
            dst = Binv[k]
            for i in range(m):
                dst[i] -= factor * src[i]
    src = Binv[row]
    for i in range(m):
        src[i] /= pivot
    basis[row] = col


def simplex(
    M: list[list[Fraction]],
    b: list[Fraction],
    N: int,
    basis: list[int],
    Binv: list[list[Fraction]],
    c: list[Fraction],
    forbidden: set[int],
    max_iter: int = 20000,
    maximize: bool = True,
) -> str:
    """Optimise ``c·z`` from a feasible basis. Mutates ``basis``/``Binv``.

    ``maximize=False`` minimises (used by phase I). Returns ``"optimal"`` or
    ``"unbounded"``. Bland's rule is used throughout to guarantee termination.
    """
    m = len(M)
    for _ in range(max_iter):
        xB = _matvec(Binv, b)
        if any(v < 0 for v in xB):
            raise SolverError("basis lost primal feasibility")
        rc = _reduced_costs(M, Binv, c, basis)
        in_basis = set(basis)
        entering = None
        for j in range(N):
            if j in forbidden or j in in_basis:
                continue
            if (rc[j] > 0) if maximize else (rc[j] < 0):
                entering = j
                break
        if entering is None:
            return "optimal"

        d = [sum((Binv[k][i] * M[i][entering] for i in range(m)), ZERO) for k in range(m)]
        best_row = None
        best_ratio = None
        for k in range(m):
            if d[k] > 0:
                ratio = xB[k] / d[k]
                if (
                    best_ratio is None
                    or ratio < best_ratio
                    or (ratio == best_ratio and basis[k] < basis[best_row])
                ):
                    best_ratio = ratio
                    best_row = k
        if best_row is None:
            return "unbounded"
        _pivot(Binv, basis, best_row, entering, d)
    raise SolverError("simplex did not terminate within iteration cap")


def _basis_values(
    M: list[list[Fraction]], b: list[Fraction], basis: list[int], Binv: list[list[Fraction]], N: int
) -> list[Fraction]:
    xB = _matvec(Binv, b)
    z = [ZERO] * N
    for k, col in enumerate(basis):
        z[col] = xB[k]
    return z


def _phase_one(
    M: list[list[Fraction]],
    b: list[Fraction],
    N: int,
    basis: list[int],
    artificial: set[int],
):
    Binv = [[Fraction(1) if i == k else ZERO for i in range(len(M))] for k in range(len(M))]
    c = [Fraction(1) if j in artificial else ZERO for j in range(N)]
    status = simplex(M, b, N, basis, Binv, c, forbidden=artificial, maximize=False)
    if status != "optimal":
        raise SolverError("phase one unexpectedly unbounded")
    z = _basis_values(M, b, basis, Binv, N)
    residual = sum((z[j] for j in artificial), ZERO)
    if residual != 0:
        return None, basis, Binv
    # Drive artificials out of the basis (all have value zero here). On a
    # degenerate row any nonzero direction entry gives a valid basis change that
    # leaves every basic value unchanged, so we do not require a positive one.
    # A row whose canonical form has no non-artificial entry at all is ``0 = 0``
    # and is dropped.
    in_basis = set(basis)
    redundant: list[int] = []
    for r in range(len(M)):
        if basis[r] not in artificial:
            continue
        row_canon = _row_canonical(Binv, M, r)
        pivot_col = None
        for j in range(N):
            if j in artificial or j in in_basis:
                continue
            if row_canon[j] != 0:
                pivot_col = j
                break
        if pivot_col is None:
            redundant.append(r)
            continue
        d = [sum((Binv[k][i] * M[i][pivot_col] for i in range(len(M))), ZERO) for k in range(len(M))]
        if d[r] == 0:  # pragma: no cover - canonical entry implies nonzero
            redundant.append(r)
            continue
        _pivot(Binv, basis, r, pivot_col, d)
        in_basis = set(basis)

    # Drop rows that reduced to ``0 = 0`` during phase I.
    if redundant:
        keep = [r for r in range(len(M)) if r not in redundant]
        M[:] = [[M[r][j] for j in range(N)] for r in keep]
        b[:] = [b[r] for r in keep]
        basis[:] = [basis[r] for r in keep]
        Binv[:] = [[Binv[r][i] for i in keep] for r in keep]
    z = _basis_values(M, b, basis, Binv, N)
    return z, basis, Binv


def _verify(
    rows: Sequence[tuple[Sequence[Fraction], str, Fraction]], x: Sequence[Fraction]
) -> None:
    for coeffs, op, r in rows:
        lhs = sum((c * v for c, v in zip(coeffs, x)), ZERO)
        ok = {
            "<=": lhs <= r,
            ">=": lhs >= r,
            "=": lhs == r,
        }[op]
        if not ok:
            raise SolverError(
                f"internal exactness check failed: {lhs} {op} {r} violated at {[qstr(v) for v in x]}"
            )


def solve(
    rows: Sequence[tuple[Sequence[Fraction], str, Fraction]],
    n: int,
    bounds: Sequence[tuple[Fraction | None, Fraction | None]],
    objectives: Iterable[int] = (),
) -> dict:
    """Feasibility plus exact min/max for each requested variable index.

    Returns a dict:
      ``feasible``   bool
      ``point``      one exact feasible vertex (or ``None``)
      ``endpoints``  ``{index: (min, max)}`` as exact ``Fraction``
      ``optima``     ``{index: (argmin_point, argmax_point)}``

    Phase I is run once; each of the ``2 * len(objectives)`` objectives then
    starts from a *copy* of that feasible basis. Endpoints are exact vertices.
    """
    objectives = list(objectives)
    M, b, basis, artificial, N = build_standard_form(rows, n, bounds)
    z, basis, Binv = _phase_one(M, b, N, list(basis), artificial)
    if z is None:
        return {"feasible": False, "point": None, "endpoints": {}, "optima": {}}
    point = z[:n]
    _verify(rows, point)

    def copy_state():
        return (
            [row[:] for row in M],
            b[:],
            basis[:],
            [row[:] for row in Binv],
        )

    endpoints: dict[int, tuple[Fraction, Fraction]] = {}
    optima: dict[int, tuple[list[Fraction], list[Fraction]]] = {}
    for idx in objectives:
        results = {}
        for sense in ("max", "min"):
            M2, b2, basis2, Binv2 = copy_state()
            # minimise c.z : max x => c = -e_idx, min x => c = +e_idx
            c = [ZERO] * N
            c[idx] = Fraction(-1) if sense == "max" else Fraction(1)
            status = simplex(M2, b2, N, basis2, Binv2, c, forbidden=artificial, maximize=False)
            if status == "unbounded":
                raise Unbounded(f"variable {idx} is unbounded")
            zopt = _basis_values(M2, b2, basis2, Binv2, N)
            xopt = zopt[:n]
            _verify(rows, xopt)
            results[sense] = (xopt[idx], xopt)
        vmin, pmin = results["min"]
        vmax, pmax = results["max"]
        endpoints[idx] = (vmin, vmax)
        optima[idx] = (pmin, pmax)
    return {"feasible": True, "point": point, "endpoints": endpoints, "optima": optima}
