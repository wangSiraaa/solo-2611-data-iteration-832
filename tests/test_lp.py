import random
from fractions import Fraction as F
from itertools import combinations

from prepquant.lp import solve

B = 5


def _det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _cramer(A, rhs):
    d = _det3(A)
    if d == 0:
        return None
    out = []
    for c in range(3):
        M = [row[:] for row in A]
        for r in range(3):
            M[r][c] = rhs[r]
        out.append(_det3(M) / d)
    return tuple(out)


def brute(rows):
    n = 3
    lines = [(list(c), rhs) for c, op, rhs in rows]
    for j in range(n):
        e = [F(0)] * n
        e[j] = F(1)
        lines.append((e[:], F(0)))
        lines.append((e[:], F(B)))
    points = set()
    for combo in combinations(range(len(lines)), n):
        p = _cramer([lines[i][0] for i in combo], [lines[i][1] for i in combo])
        if p and all(0 <= v <= B for v in p):
            points.add(p)
    feasible = []
    for x in points:
        if all(
            {
                "<=": sum(c[i] * x[i] for i in range(n)) <= rhs,
                ">=": sum(c[i] * x[i] for i in range(n)) >= rhs,
                "=": sum(c[i] * x[i] for i in range(n)) == rhs,
            }[op]
            for c, op, rhs in rows
        ):
            feasible.append(x)
    return feasible


def test_lp_matches_bruteforce():
    rng = random.Random(2026)
    for _ in range(150):
        rows = []
        for _ in range(rng.randint(1, 6)):
            c = [F(rng.randint(-2, 2), rng.choice([1, 2])) for _ in range(3)]
            if all(v == 0 for v in c):
                c[0] = F(1)
            rhs = F(rng.randint(-3, 6), rng.choice([1, 2]))
            rows.append((c, rng.choice(["<=", ">=", "="]), rhs))
        bounds = [(F(0), F(B))] * 3
        got = solve(rows, 3, bounds, [0, 1, 2])
        expected = brute(rows)
        assert got["feasible"] == bool(expected), (rows, expected, got)
        if expected:
            for i in range(3):
                lo = min(p[i] for p in expected)
                hi = max(p[i] for p in expected)
                assert got["endpoints"][i] == (lo, hi), (i, rows, expected, got)


def test_zero_variable_feasibility():
    assert solve([([], ">=", 5)], 0, [], [])["feasible"] is False
    assert solve([([], "<=", 4)], 0, [], [])["feasible"] is True


def test_lower_bounds():
    r = solve([], 1, [(F(3), F(7))], [0])
    assert r["endpoints"][0] == (F(3), F(7))
    r = solve([([F(1), F(1)], "<=", F(10))], 2, [(F(2), F(100)), (F(5), F(100))], [0, 1])
    assert r["endpoints"][0] == (F(2), F(5))
    assert r["endpoints"][1] == (F(5), F(8))
