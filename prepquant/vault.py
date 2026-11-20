"""Quantitation vault: measurements, linear interval constraints, exact inversion.

The vault is deliberately independent of the ledger's private tables. It accepts
*product payloads* in the exported format (volume plus a mass expression in the
sample unknowns) and stores a snapshot of each payload, so an analysis can always
be reproduced even after the ledger moves on.

A measurement is a concentration reading of one product:
``observed ± uncertainty`` (mg/mL). The truth is constrained to
``[observed - uncertainty, observed + uncertainty] ∩ [0, ∞)``. Since
``concentration = mass_expression / volume`` and ``volume > 0``, this becomes the
two exact linear rows

    mass_expression >= max(0, observed - uncertainty) * volume
    mass_expression <= (observed + uncertainty) * volume

Every analysis solves *all* its measurements jointly over the given per-sample
bounds. Endpoints are produced by the exact simplex in :mod:`prepquant.lp`; no
averaging is ever applied to a non-unique solution.
"""

from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .expr import LinExpr
from .ledger import MAX_SAMPLES, LedgerError, Product, product_from_payload
from .lp import Infeasible, SolverError, solve as lp_solve
from .rational import ZERO, nonneg, positive, q, qstr

MAX_MEASUREMENTS = 20


class VaultError(Exception):
    pass


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class Vault:
    def __init__(self, db_path: str | Path) -> None:
        self.db_path = str(db_path)
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._create_schema()

    def _create_schema(self) -> None:
        with self.conn:
            self.conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS measurements(
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL UNIQUE,
                    product_id TEXT NOT NULL,
                    payload TEXT NOT NULL,
                    concentration TEXT NOT NULL,
                    observed TEXT NOT NULL,
                    uncertainty TEXT NOT NULL,
                    created_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS analyses(
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL,
                    status TEXT NOT NULL,
                    measurement_ids TEXT NOT NULL,
                    bounds TEXT NOT NULL,
                    result TEXT NOT NULL,
                    created_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS analysis_measurements(
                    analysis_id INTEGER NOT NULL,
                    measurement_id INTEGER NOT NULL,
                    product_id TEXT NOT NULL,
                    PRIMARY KEY(analysis_id, measurement_id)
                );
                """
            )

    def close(self) -> None:
        self.conn.close()

    # ------------------------------------------------------------ measurements
    def add_measurement(
        self,
        name: str,
        product: Mapping | Product,
        observed,
        uncertainty,
    ) -> dict:
        if not isinstance(name, str) or not name:
            raise VaultError("measurement name must be a non-empty string")
        try:
            payload = product.to_payload() if isinstance(product, Product) else dict(product)
        except AttributeError as exc:
            raise VaultError(f"unusable product: {exc}") from exc

        try:
            prod = product_from_payload(payload)
        except LedgerError as exc:
            raise VaultError(str(exc)) from exc
        if prod.volume <= 0:
            raise VaultError(f"product {prod.product_id} has non-positive volume {qstr(prod.volume)}")
        if len(prod.unknowns) > MAX_SAMPLES:
            raise VaultError(f"product references {len(prod.unknowns)} unknowns (max {MAX_SAMPLES})")

        obs = nonneg(q(observed, "observed"), "observed")
        unc = nonneg(q(uncertainty, "uncertainty"), "uncertainty")
        conc = prod.mass.scale(Fraction(1) / prod.volume)

        try:
            with self.conn:
                cur = self.conn.execute(
                    "INSERT INTO measurements(name, product_id, payload, concentration, observed, uncertainty, created_at)"
                    " VALUES (?,?,?,?,?,?,?)",
                    (name, prod.product_id, json.dumps(payload, sort_keys=True), conc.formatted(),
                     qstr(obs), qstr(unc), _now()),
                )
        except sqlite3.IntegrityError as exc:
            raise VaultError(f"measurement {name!r} already exists; measurements are immutable") from exc
        return self.measurement(int(cur.lastrowid))

    def measurement(self, ref: str | int) -> dict:
        row = self._measurement_row(ref)
        return self._measurement_dict(row)

    def _measurement_row(self, ref: str | int) -> sqlite3.Row:
        if isinstance(ref, int):
            row = self.conn.execute("SELECT * FROM measurements WHERE id=?", (ref,)).fetchone()
        else:
            row = self.conn.execute("SELECT * FROM measurements WHERE name=?", (str(ref),)).fetchone()
        if row is None:
            raise VaultError(f"measurement {ref!r} not found")
        return row

    def _measurement_dict(self, row: sqlite3.Row) -> dict:
        payload = json.loads(row["payload"])
        prod = product_from_payload(payload)
        obs = q(row["observed"], "observed")
        unc = q(row["uncertainty"], "uncertainty")
        lo = max(ZERO, obs - unc)
        hi = obs + unc
        return {
            "id": row["id"],
            "name": row["name"],
            "product_id": row["product_id"],
            "conc_expr": prod.mass.scale(Fraction(1) / prod.volume),
            "mass_expr": prod.mass,
            "volume": prod.volume,
            "observed": obs,
            "uncertainty": unc,
            "interval": (lo, hi),
            "payload": payload,
            "unknowns": prod.unknowns,
            "created_at": row["created_at"],
        }

    def list_measurements(self) -> list[dict]:
        rows = self.conn.execute("SELECT * FROM measurements ORDER BY id").fetchall()
        out = []
        for row in rows:
            m = self._measurement_dict(row)
            m["conc_expr"] = m["conc_expr"].formatted()
            m["observed"] = qstr(m["observed"])
            m["uncertainty"] = qstr(m["uncertainty"])
            m["interval"] = [qstr(m["interval"][0]), qstr(m["interval"][1])]
            m.pop("payload", None)
            out.append(m)
        return out

    # ---------------------------------------------------------------- analysis
    def analyze(
        self,
        name: str,
        measurement_refs: Sequence[str | int],
        bounds: Mapping[str, tuple],
        ledger=None,
        record_consumption: bool = True,
    ) -> dict:
        if not isinstance(name, str) or not name:
            raise VaultError("analysis name must be a non-empty string")
        refs = list(measurement_refs)
        if not refs:
            raise VaultError("an analysis needs at least one measurement")
        if len(refs) > MAX_MEASUREMENTS:
            raise VaultError(f"at most {MAX_MEASUREMENTS} measurements per analysis")
        seen_refs: set[int] = set()
        measurements = []
        for ref in refs:
            m = self.measurement(ref)
            if m["id"] in seen_refs:
                raise VaultError(f"measurement {m['name']!r} listed more than once")
            seen_refs.add(m["id"])
            measurements.append(m)

        # Global unknown ordering: first-seen across measurement products.
        unknowns: list[str] = []
        for m in measurements:
            for u in m["unknowns"]:
                if u not in unknowns:
                    unknowns.append(u)
        if len(unknowns) > MAX_SAMPLES:
            raise VaultError(f"{len(unknowns)} sample unknowns exceed the limit of {MAX_SAMPLES}")

        parsed_bounds: dict[str, tuple[Fraction, Fraction]] = {}
        for u in unknowns:
            if u not in bounds:
                raise VaultError(f"missing bounds for sample {u!r}")
            lo, hi = bounds[u]
            lo = nonneg(q(lo, f"bound[{u}].lo"), f"bound[{u}].lo")
            hi = nonneg(q(hi, f"bound[{u}].hi"), f"bound[{u}].hi")
            if lo > hi:
                raise VaultError(f"lower bound {qstr(lo)} exceeds upper bound {qstr(hi)} for {u!r}")
            parsed_bounds[u] = (lo, hi)
        # Bounds for samples that no measured product references are harmless:
        # they do not constrain this analysis and are not solved for.

        rows: list[tuple[list[Fraction], str, Fraction]] = []
        interval_rows = []
        for m in measurements:
            # Constraints are written on *mass* (mg), which stays linear no matter
            # how many top-ups and mixes the product went through:
            #   lo*V <= mass_const + Σ mass_coeff*x <= hi*V
            mass = m["mass_expr"]
            coeffs = [mass.coeff_of(u) for u in unknowns]
            const = mass.const
            lo, hi = m["interval"]
            rows.append((coeffs, "<=", hi * m["volume"] - const))
            rows.append((coeffs, ">=", lo * m["volume"] - const))
            interval_rows.append(
                {
                    "measurement": m["name"],
                    "product_id": m["product_id"],
                    "concentration": m["conc_expr"].formatted(),
                    "mass": mass.formatted(),
                    "volume": qstr(m["volume"]),
                    "interval": [qstr(lo), qstr(hi)],
                    "mass_coeffs": {u: qstr(c) for u, c in zip(unknowns, coeffs)},
                    "mass_constant": qstr(const),
                }
            )

        index = {u: i for i, u in enumerate(unknowns)}
        solved = lp_solve(
            rows,
            len(unknowns),
            [parsed_bounds[u] for u in unknowns],
            [index[u] for u in unknowns],
        )

        individual_conflicts = self._individual_conflicts(measurements, unknowns, parsed_bounds)

        if not solved["feasible"]:
            status = "contradiction"
            result = {
                "status": status,
                "unknowns": unknowns,
                "bounds": {u: [qstr(lo), qstr(hi)] for u, (lo, hi) in parsed_bounds.items()},
                "constraints": interval_rows,
                "measurements": [m["name"] for m in measurements],
                "individual_conflicts": individual_conflicts,
            }
        else:
            endpoints = {}
            unique_all = True
            point = solved["point"]
            for u in unknowns:
                lo, hi = solved["endpoints"][index[u]]
                unique = lo == hi
                unique_all = unique_all and unique
                endpoints[u] = {"min": qstr(lo), "max": qstr(hi), "unique": unique}
            status = "feasible"
            result = {
                "status": status,
                "unknowns": unknowns,
                "bounds": {u: [qstr(lo), qstr(hi)] for u, (lo, hi) in parsed_bounds.items()},
                "constraints": interval_rows,
                "measurements": [m["name"] for m in measurements],
                "witness_point": {u: qstr(point[index[u]]) for u in unknowns},
                "endpoints": endpoints,
                "endpoint_points": {
                    u: {
                        "min": {v: qstr(solved["optima"][index[u]][0][index[v]]) for v in unknowns},
                        "max": {v: qstr(solved["optima"][index[u]][1][index[v]]) for v in unknowns},
                    }
                    for u in unknowns
                },
                "unique": unique_all,
                "individual_conflicts": individual_conflicts,
            }

        # Decide consumption before writing: ids absent from the ledger are
        # external products with no version to protect, so they are simply not
        # consumed. This keeps a mixed or purely external analysis from failing
        # after the analysis row has been written.
        to_consume: list[str] = []
        if ledger is not None and record_consumption:
            for m in measurements:
                if ledger.has_product(m["product_id"]):
                    to_consume.append(m["product_id"])

        # persist the analysis and its expression snapshot (old rows are never touched)
        snapshot = []
        for m in measurements:
            snapshot.append(
                {
                    "name": m["name"],
                    "product_id": m["product_id"],
                    "concentration": str(m["conc_expr"]),
                    "interval": [qstr(m["interval"][0]), qstr(m["interval"][1])],
                    "volume": qstr(m["volume"]),
                    "payload": m["payload"],
                }
            )
        result["measurement_snapshot"] = snapshot
        result["name"] = name
        result["created_at"] = _now()
        with self.conn:
            cur = self.conn.execute(
                "INSERT INTO analyses(name, status, measurement_ids, bounds, result, created_at) VALUES (?,?,?,?,?,?)",
                (
                    name,
                    status,
                    json.dumps([m["id"] for m in measurements]),
                    json.dumps({u: [qstr(lo), qstr(hi)] for u, (lo, hi) in parsed_bounds.items()}, sort_keys=True),
                    json.dumps(result, sort_keys=True),
                    result["created_at"],
                ),
            )
            aid = int(cur.lastrowid)
            for m in measurements:
                self.conn.execute(
                    "INSERT INTO analysis_measurements(analysis_id, measurement_id, product_id) VALUES (?,?,?)",
                    (aid, m["id"], m["product_id"]),
                )
        result["id"] = aid

        # Consume the pinned product versions in the ledger. Consumption is
        # recorded only for analyses that actually reference the payload; old
        # analyses and their pinned snapshots are untouched.
        for pid in to_consume:
            ledger.mark_consumed(pid, ref=f"analysis:{aid}:{name}")
        return result

    def _individual_conflicts(self, measurements, unknowns, bounds) -> list[dict]:
        """A measurement alone is contradictory if its interval misses the mass
        range attainable over the sample bounds. This is exact and cheap."""
        conflicts = []
        for m in measurements:
            expr = m["mass_expr"]
            mass_min = expr.const
            mass_max = expr.const
            for u in unknowns:
                c = expr.coeff_of(u)
                lo, hi = bounds[u]
                if c >= 0:
                    mass_min += c * lo
                    mass_max += c * hi
                else:
                    mass_min += c * hi
                    mass_max += c * lo
            lo_c, hi_c = m["interval"]
            if hi_c * m["volume"] < mass_min or lo_c * m["volume"] > mass_max:
                conflicts.append(
                    {
                        "measurement": m["name"],
                        "reason": "reading interval is disjoint from the mass attainable within bounds",
                        "attainable_mass": [qstr(mass_min), qstr(mass_max)],
                        "reading_mass_interval": [qstr(lo_c * m["volume"]), qstr(hi_c * m["volume"])],
                    }
                )
        return conflicts

    # -------------------------------------------------------------- read side
    def list_analyses(self) -> list[dict]:
        rows = self.conn.execute("SELECT id, name, status, measurement_ids, created_at FROM analyses ORDER BY id").fetchall()
        return [
            {
                "id": r["id"],
                "name": r["name"],
                "status": r["status"],
                "measurement_ids": json.loads(r["measurement_ids"]),
                "created_at": r["created_at"],
            }
            for r in rows
        ]

    def analysis(self, ref: int | str) -> dict:
        if isinstance(ref, int):
            row = self.conn.execute("SELECT * FROM analyses WHERE id=?", (ref,)).fetchone()
        else:
            row = self.conn.execute("SELECT * FROM analyses WHERE name=? ORDER BY id DESC", (str(ref),)).fetchone()
        if row is None:
            raise VaultError(f"analysis {ref!r} not found")
        result = json.loads(row["result"])
        result["id"] = row["id"]
        result["name"] = row["name"]
        return result

    def consumed_product_ids(self) -> set[str]:
        rows = self.conn.execute("SELECT DISTINCT product_id FROM analysis_measurements").fetchall()
        return {r["product_id"] for r in rows}
