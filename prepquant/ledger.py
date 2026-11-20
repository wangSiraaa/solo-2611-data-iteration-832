"""Preparation ledger: an ordered, atomic, replayable log of lab operations.

Design choices and their consequences
-------------------------------------
* **Replay, not mutation.** The durable truth is the ordered operation log.
  Container state, product versions and provenance are rebuilt from that log.
  A correction rewrites the log from the state *before* the corrected step and
  reproduces every later step; if the replay fails, nothing is committed. This
  makes the "whole batch is atomic" and "correction replays from a prior state"
  guarantees structural rather than best-effort.
* **Content-addressed product versions.** A product id is
  ``P<seq>-<container>-<content hash>``. Analyses can only pin a product id, so
  a correction that would silently change the contents behind a consumed id is
  refused, while an unconsumed id simply disappears and a stale export cannot be
  imported by accident.
* **Mass is the currency.** Every container stores ``mass = const + Σ a_i x_i``
  (mg) and a constant volume (mL). Concentration is ``mass / volume``. Unknowns
  are the original sample concentrations. All arithmetic is :class:`Fraction`.
* **Provenance is first-class.** Each container carries an *origin book* mapping
  the originating container (a standard batch, a sample batch, or ``solvent``)
  to the volume and mass expression that origin still contributes. Transfer
  splits an origin pro-rata by volume; mix unions books; solvent adds a zero-mass
  origin; evaporation rescales volume only (solute mass is untouched); discard
  empties the book. No measurement is needed to obtain this.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from dataclasses import dataclass, field
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .expr import LinExpr
from .rational import ZERO, nonneg, positive, q, qstr

MAX_SAMPLES = 6
MAX_CONTAINERS = 80
SOLVENT_ORIGIN = "solvent"

OP_KINDS = ("standard", "sample", "transfer", "mix", "solvent", "evaporate", "discard")


class LedgerError(Exception):
    """A batch violated the physical or bookkeeping rules; nothing changed."""


@dataclass
class Container:
    cid: str
    kind: str
    volume: Fraction
    mass: LinExpr
    origins: dict[str, tuple[Fraction, LinExpr]] = field(default_factory=dict)
    sample: str | None = None
    label: str | None = None
    version: int = 0

    def clone(self) -> "Container":
        return Container(
            cid=self.cid,
            kind=self.kind,
            volume=self.volume,
            mass=LinExpr(self.mass.const, self.mass.coeffs, self.mass.unknowns),
            origins={k: (v, LinExpr(m.const, m.coeffs, m.unknowns)) for k, (v, m) in self.origins.items()},
            sample=self.sample,
            label=self.label,
            version=self.version,
        )

    def is_empty(self) -> bool:
        return self.volume == 0 and self.mass.is_constant() and self.mass.const == 0

    def normalize_origins(self) -> None:
        self.origins = {
            k: (vol, mass)
            for k, (vol, mass) in self.origins.items()
            if vol != 0 or not mass.is_constant() or mass.const != 0
        }


def _zero_expr(unknowns: Sequence[str]) -> LinExpr:
    return LinExpr.constant(ZERO, unknowns)


def _check_keys(payload: Mapping, required: Iterable[str], optional: Iterable[str] = (), what: str = "operation") -> None:
    required = set(required)
    optional = set(optional)
    missing = required - set(payload)
    if missing:
        raise LedgerError(f"{what} is missing field(s): {', '.join(sorted(missing))}")
    extra = set(payload) - required - optional
    if extra:
        raise LedgerError(f"{what} has unknown field(s): {', '.join(sorted(extra))}")


def op_kind(op: Mapping) -> str:
    kind = op.get("op")
    if kind not in OP_KINDS:
        raise LedgerError(f"unknown operation {kind!r}; expected one of {', '.join(OP_KINDS)}")
    return kind


def sample_names_in_order(ops: Sequence[Mapping]) -> list[str]:
    seen: list[str] = []
    for op in ops:
        kind = op_kind(op)
        if kind == "sample":
            name = op.get("sample")
            if not isinstance(name, str) or not name:
                raise LedgerError("sample operation requires a non-empty 'sample' name")
            if name not in seen:
                seen.append(name)
    if len(seen) > MAX_SAMPLES:
        raise LedgerError(f"at most {MAX_SAMPLES} sample unknowns are supported, got {len(seen)}")
    return seen


def _referenced_unknowns(container: Container) -> list[str]:
    """Samples that actually appear with a non-zero coefficient.

    Trimming keeps exported products minimal and, crucially, keeps product ids
    stable when a *later* step declares an additional sample: no zero coefficient
    ever enters the identity of an existing product.
    """
    used: list[str] = []

    def scan(expr: LinExpr) -> None:
        for name, coeff in expr.sparse().items():
            if name not in used:
                used.append(name)

    scan(container.mass)
    for _, mass in container.origins.values():
        scan(mass)
    # keep the ledger's declaration order rather than discovery order
    order = {name: i for i, name in enumerate(container.mass.unknowns)}
    return sorted(used, key=lambda n: order.get(n, 10**9))


def _expr_terms(expr: LinExpr) -> dict:
    return {"const": qstr(expr.const), "terms": {k: qstr(v) for k, v in expr.sparse().items()}}


def _product_fingerprint(container: Container) -> str:
    payload = {
        "cid": container.cid,
        "kind": container.kind,
        "volume": qstr(container.volume),
        "mass": _expr_terms(container.mass),
        "origins": {
            k: {"volume": qstr(v), "mass": _expr_terms(m)}
            for k, (v, m) in sorted(container.origins.items())
        },
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()[:16]


@dataclass
class Product:
    product_id: str
    container_id: str
    version: int
    seq: int
    kind: str
    volume: Fraction
    mass: LinExpr
    origins: dict[str, tuple[Fraction, LinExpr]]
    unknowns: tuple[str, ...]

    def referenced_unknowns(self) -> list[str]:
        used: list[str] = []
        for expr in [self.mass] + [m for _, m in self.origins.values()]:
            for name in expr.sparse():
                if name not in used:
                    used.append(name)
        order = {name: i for i, name in enumerate(self.unknowns)}
        return sorted(used, key=lambda n: order.get(n, 10**9))

    def to_payload(self) -> dict:
        """External format: volume plus a mass expression over *referenced* samples.

        Zero-coefficient samples are omitted so the payload can be solved without
        inventing bounds for samples the product does not actually contain.
        """
        return {
            "product_id": self.product_id,
            "container": self.container_id,
            "version": self.version,
            "seq": self.seq,
            "kind": self.kind,
            "volume": qstr(self.volume),
            "unknowns": self.referenced_unknowns(),
            "mass": _expr_terms(self.mass),
            "origins": [
                {
                    "origin": origin,
                    "volume": qstr(vol),
                    "mass": _expr_terms(mass),
                }
                for origin, (vol, mass) in sorted(self.origins.items())
            ],
        }


class Ledger:
    def __init__(self, db_path: str | Path) -> None:
        self.db_path = str(db_path)
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._create_schema()
        self.ops: list[dict] = []
        self.containers: dict[str, Container] = {}
        self.products: list[Product] = []
        self.unknowns: tuple[str, ...] = ()
        self._reload()

    # ------------------------------------------------------------------ db
    def _create_schema(self) -> None:
        with self.conn:
            self.conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS ops(seq INTEGER PRIMARY KEY, kind TEXT NOT NULL, payload TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS corrections(
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    target_seq INTEGER NOT NULL,
                    replaced_payload TEXT NOT NULL,
                    new_payload TEXT NOT NULL,
                    note TEXT,
                    created_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS products(
                    product_id TEXT PRIMARY KEY,
                    container_id TEXT NOT NULL,
                    version INTEGER NOT NULL,
                    seq INTEGER NOT NULL,
                    kind TEXT NOT NULL,
                    volume TEXT NOT NULL,
                    mass TEXT NOT NULL,
                    origins TEXT NOT NULL,
                    unknowns TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS consumed(
                    product_id TEXT NOT NULL,
                    ref TEXT NOT NULL,
                    consumed_at TEXT NOT NULL,
                    PRIMARY KEY(product_id, ref)
                );
                """
            )

    def close(self) -> None:
        self.conn.close()

    # --------------------------------------------------------------- replay
    def _reload(self) -> None:
        rows = self.conn.execute("SELECT seq, kind, payload FROM ops ORDER BY seq").fetchall()
        self.ops = [{"seq": r["seq"], "kind": r["kind"], "op": json.loads(r["payload"])} for r in rows]
        self._replay(self.ops, persist=False)

    def _replay(self, ops: Sequence[Mapping], persist: bool) -> None:
        """Replay ``ops`` from an empty state.

        ``ops`` items are ``{"seq", "kind", "op"}``. Raises :class:`LedgerError`
        on any violation, leaving the caller's stored state untouched.
        """
        payloads = [record["op"] for record in ops]
        for record in ops:
            if op_kind(record["op"]) != record["kind"]:
                raise LedgerError(f"seq {record['seq']}: stored kind does not match payload")
        try:
            unknowns = tuple(sample_names_in_order(payloads))
        except LedgerError:
            raise
        except (ValueError, TypeError) as exc:
            raise LedgerError(str(exc)) from exc
        containers: dict[str, Container] = {}
        products: list[Product] = []
        self._replay_products = products

        for record in ops:
            seq = record["seq"]
            op = record["op"]
            kind = record["kind"]
            try:
                self._replay_one(seq, kind, op, containers, unknowns)
            except LedgerError:
                raise
            except (ValueError, TypeError) as exc:
                raise LedgerError(f"seq {seq}: {exc}") from exc
            if len(containers) > MAX_CONTAINERS:
                raise LedgerError(f"at most {MAX_CONTAINERS} containers are supported")

        self.containers = containers
        self.unknowns = unknowns
        self.products = products
        if persist:
            self._persist_products(products)


    def _replay_one(self, seq, kind, op, containers, unknowns) -> None:
        handler = {
            "standard": self._op_standard,
            "sample": self._op_sample,
            "transfer": self._op_transfer,
            "mix": self._op_mix,
            "solvent": self._op_solvent,
            "evaporate": self._op_evaporate,
            "discard": self._op_discard,
        }.get(kind)
        if handler is None:  # pragma: no cover - guarded earlier
            raise LedgerError(f"unknown operation {kind!r}")
        handler(seq, op, containers, unknowns)

    def _validate_common(self, seq: int, containers: Mapping[str, Container], cid: str) -> Container:
        if cid not in containers:
            raise LedgerError(f"seq {seq}: container {cid!r} does not exist")
        return containers[cid]

    def _check_new_container_id(self, seq: int, containers: Mapping[str, Container], cid: str) -> None:
        if not isinstance(cid, str) or not cid:
            raise LedgerError(f"seq {seq}: container id must be a non-empty string")
        if cid in containers and not containers[cid].is_empty():
            raise LedgerError(f"seq {seq}: container {cid!r} already holds liquid")

    def _record(self, seq: int, kind: str, container: Container, unknowns: Sequence[str]) -> None:
        container.version += 1
        container.normalize_origins()
        product_id = f"P{seq}-{container.cid}-{_product_fingerprint(container)}"
        self._replay_products.append(
            Product(
                product_id=product_id,
                container_id=container.cid,
                version=container.version,
                seq=seq,
                kind=container.kind,
                volume=container.volume,
                mass=container.mass,
                origins=dict(container.origins),
                unknowns=tuple(unknowns),
            )
        )

    def _op_standard(self, seq, op, containers, unknowns) -> None:
        _check_keys(op, ("op", "container", "volume", "concentration"), ("label",), "standard")
        cid = op["container"]
        self._check_new_container_id(seq, containers, cid)
        volume = positive(q(op["volume"], "volume"), "volume")
        conc = nonneg(q(op["concentration"], "concentration"), "concentration")
        version = 0
        if cid in containers:
            if not containers[cid].is_empty():
                raise LedgerError(f"seq {seq}: container {cid!r} already holds liquid; container identity is stable")
            version = containers[cid].version
        mass = LinExpr.constant(volume * conc, unknowns)
        book = {cid: (volume, mass)}
        container = Container(cid, "standard", volume, mass, book, label=op.get("label"), version=version)
        containers[cid] = container
        self._record(seq, "standard", container, unknowns)

    def _op_sample(self, seq, op, containers, unknowns) -> None:
        _check_keys(op, ("op", "container", "volume", "sample"), ("label",), "sample")
        cid = op["container"]
        self._check_new_container_id(seq, containers, cid)
        name = op["sample"]
        if name not in unknowns:
            raise LedgerError(f"seq {seq}: sample {name!r} not declared in the unknown ordering")
        volume = positive(q(op["volume"], "volume"), "volume")
        version = 0
        if cid in containers:
            if not containers[cid].is_empty():
                raise LedgerError(f"seq {seq}: container {cid!r} already holds liquid; container identity is stable")
            version = containers[cid].version
        mass = LinExpr.variable(name, unknowns, scale=volume)
        book = {cid: (volume, mass)}
        container = Container(cid, "sample", volume, mass, book, sample=name, label=op.get("label"), version=version)
        containers[cid] = container
        self._record(seq, "sample", container, unknowns)

    def _op_transfer(self, seq, op, containers, unknowns) -> None:
        _check_keys(op, ("op", "source", "volume", "destination"), ("label",), "transfer")
        src = self._validate_common(seq, containers, op["source"])
        dst_id = op["destination"]
        if not isinstance(dst_id, str) or not dst_id:
            raise LedgerError(f"seq {seq}: destination must be a non-empty string")
        if dst_id == src.cid:
            raise LedgerError(f"seq {seq}: cannot transfer into the source container itself")
        volume = positive(q(op["volume"], "volume"), "volume")
        if volume > src.volume:
            raise LedgerError(
                f"seq {seq}: cannot take {qstr(volume)} mL from {src.cid} which holds only {qstr(src.volume)} mL"
            )
        if dst_id in containers and not containers[dst_id].is_empty():
            raise LedgerError(f"seq {seq}: destination {dst_id!r} is not empty")
        frac = volume / src.volume
        moved: dict[str, tuple[Fraction, LinExpr]] = {}
        for origin, (ov, om) in src.origins.items():
            mv = ov * frac
            mm = om.scale(frac)
            moved[origin] = (mv, mm)
        # shrink source by the exact remainder factor
        remain = (src.volume - volume) / src.volume
        src.origins = {o: (ov * remain, om.scale(remain)) for o, (ov, om) in src.origins.items()}
        src.volume -= volume
        src.mass = LinExpr.constant(ZERO, unknowns)
        for _, (_, m) in src.origins.items():
            src.mass = src.mass + m
        src.normalize_origins()
        self._record(seq, "transfer", src, unknowns)

        if dst_id in containers:
            dst = containers[dst_id]
            dst.kind = "mixture"
        else:
            dst = Container(dst_id, "mixture", ZERO, LinExpr.constant(ZERO, unknowns))
            containers[dst_id] = dst
        dst.volume = volume
        dst.mass = LinExpr.constant(ZERO, unknowns)
        for _, (_, m) in moved.items():
            dst.mass = dst.mass + m
        dst.origins = dict(moved)
        self._record(seq, "transfer", dst, unknowns)

    def _op_mix(self, seq, op, containers, unknowns) -> None:
        _check_keys(op, ("op", "sources", "destination"), ("label",), "mix")
        sources = op["sources"]
        if not isinstance(sources, list) or not sources:
            raise LedgerError(f"seq {seq}: mix needs a non-empty 'sources' list")
        dst_id = op["destination"]
        if not isinstance(dst_id, str) or not dst_id:
            raise LedgerError(f"seq {seq}: destination must be a non-empty string")
        if dst_id in containers and not containers[dst_id].is_empty():
            raise LedgerError(f"seq {seq}: destination {dst_id!r} is not empty")
        requests: list[tuple[str, Fraction]] = []
        seen: set[str] = set()
        for item in sources:
            if not isinstance(item, Mapping):
                raise LedgerError(f"seq {seq}: each mix source must be an object")
            _check_keys(item, ("container", "volume"), what="mix source")
            cid = item["container"]
            if cid == dst_id:
                raise LedgerError(f"seq {seq}: destination cannot also be a mix source")
            if cid in seen:
                raise LedgerError(f"seq {seq}: container {cid!r} listed twice; pipe a single volume instead")
            seen.add(cid)
            volume = positive(q(item["volume"], "volume"), "volume")
            requests.append((cid, volume))

        moved_books: dict[str, tuple[Fraction, LinExpr]] = {}
        total = ZERO
        for cid, volume in requests:
            src = self._validate_common(seq, containers, cid)
            if volume > src.volume:
                raise LedgerError(
                    f"seq {seq}: cannot take {qstr(volume)} mL from {cid} which holds only {qstr(src.volume)} mL"
                )
            total += volume
            frac = volume / src.volume
            for origin, (ov, om) in src.origins.items():
                pv, pm = moved_books.get(origin, (ZERO, LinExpr.constant(ZERO, unknowns)))
                moved_books[origin] = (pv + ov * frac, pm + om.scale(frac))
            remain = (src.volume - volume) / src.volume
            src.origins = {o: (ov * remain, om.scale(remain)) for o, (ov, om) in src.origins.items()}
            src.volume -= volume
            src.mass = LinExpr.constant(ZERO, unknowns)
            for _, (_, m) in src.origins.items():
                src.mass = src.mass + m
            src.normalize_origins()
            self._record(seq, "mix", src, unknowns)

        if dst_id in containers:
            dst = containers[dst_id]
            dst.kind = "mixture"
        else:
            dst = Container(dst_id, "mixture", ZERO, LinExpr.constant(ZERO, unknowns))
            containers[dst_id] = dst
        dst.volume = total
        dst.mass = LinExpr.constant(ZERO, unknowns)
        for _, (_, m) in moved_books.items():
            dst.mass = dst.mass + m
        dst.origins = moved_books
        self._record(seq, "mix", dst, unknowns)

    def _op_solvent(self, seq, op, containers, unknowns) -> None:
        _check_keys(op, ("op", "container", "volume"), (), "solvent")
        cid = op["container"]
        if not isinstance(cid, str) or not cid:
            raise LedgerError(f"seq {seq}: container id must be a non-empty string")
        volume = positive(q(op["volume"], "volume"), "volume")
        if cid not in containers:
            # Solvent can be added to a fresh, still-empty destination first.
            containers[cid] = Container(cid, "empty", ZERO, LinExpr.constant(ZERO, unknowns))
        c = containers[cid]
        c.volume += volume
        prev_v, prev_m = c.origins.get(SOLVENT_ORIGIN, (ZERO, LinExpr.constant(ZERO, unknowns)))
        c.origins[SOLVENT_ORIGIN] = (prev_v + volume, prev_m)
        if c.kind == "empty":
            c.kind = "solvent"
        self._record(seq, "solvent", c, unknowns)

    def _op_evaporate(self, seq, op, containers, unknowns) -> None:
        _check_keys(op, ("op", "container"), ("to", "remove"), "evaporate")
        c = self._validate_common(seq, containers, op["container"])
        if ("to" in op) == ("remove" in op):
            raise LedgerError(f"seq {seq}: evaporate needs exactly one of 'to' or 'remove'")
        if c.volume <= 0:
            raise LedgerError(f"seq {seq}: cannot evaporate from empty container {c.cid}")
        if "to" in op:
            new_volume = q(op["to"], "to")
        else:
            new_volume = c.volume - q(op["remove"], "remove")
        if new_volume <= 0:
            raise LedgerError(f"seq {seq}: evaporation must leave a positive volume, got {qstr(new_volume)}")
        if new_volume >= c.volume:
            raise LedgerError(f"seq {seq}: evaporation must reduce volume ({qstr(c.volume)} -> {qstr(new_volume)})")
        factor = new_volume / c.volume
        c.origins = {o: (v * factor, m) for o, (v, m) in c.origins.items()}
        c.volume = new_volume
        # mass untouched: only solvent leaves
        self._record(seq, "evaporate", c, unknowns)

    def _op_discard(self, seq, op, containers, unknowns) -> None:
        _check_keys(op, ("op", "container"), ("reason",), "discard")
        c = self._validate_common(seq, containers, op["container"])
        if c.volume <= 0:
            raise LedgerError(f"seq {seq}: container {c.cid} is already empty")
        c.volume = ZERO
        c.mass = LinExpr.constant(ZERO, unknowns)
        c.origins = {}
        c.kind = "empty"
        self._record(seq, "discard", c, unknowns)

    # ------------------------------------------------------------- persistence
    def _persist_products(self, products: Sequence[Product]) -> None:
        with self.conn:
            self.conn.execute("DELETE FROM products")
            for p in products:
                self.conn.executemany(
                    "INSERT INTO products(product_id, container_id, version, seq, kind, volume, mass, origins, unknowns)"
                    " VALUES (?,?,?,?,?,?,?,?,?)",
                    [
                        (
                            p.product_id,
                            p.container_id,
                            p.version,
                            p.seq,
                            p.kind,
                            qstr(p.volume),
                            json.dumps(
                                {
                                    "const": qstr(p.mass.const),
                                    "terms": {k: qstr(v) for k, v in p.mass.sparse().items()},
                                }
                            ),
                            json.dumps(
                                {
                                    k: {
                                        "volume": qstr(v),
                                        "mass": {
                                            "const": qstr(m.const),
                                            "terms": {k2: qstr(v2) for k2, v2 in m.sparse().items()},
                                        },
                                    }
                                    for k, (v, m) in p.origins.items()
                                }
                            ),
                            json.dumps(list(p.unknowns)),
                        )
                    ],
                )
            self.conn.execute(
                "INSERT INTO meta(key, value) VALUES ('latest_seq', ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (str(len(self.ops)),),
            )

    # --------------------------------------------------------------- batch API
    def apply_batch(self, ops: Sequence[Mapping]) -> list[int]:
        """Validate and append an ordered batch atomically. Returns new seq ids."""
        if not ops:
            return []
        candidate = [dict(o) for o in self.ops]
        next_seq = len(candidate) + 1
        for op in ops:
            payload = json.loads(json.dumps(op))  # deep copy, JSON-normalised
            candidate.append({"seq": next_seq, "kind": op_kind(payload), "op": payload})
            next_seq += 1
        # replay validates every rule; any exception means nothing is written
        self._replay(candidate, persist=False)
        # commit
        with self.conn:
            for record in candidate[len(self.ops):]:
                self.conn.execute(
                    "INSERT INTO ops(seq, kind, payload) VALUES (?,?,?)",
                    (record["seq"], record["kind"], json.dumps(record["op"], sort_keys=True)),
                )
        self.ops = candidate
        self._replay(candidate, persist=False)
        self._persist_products(self.products)
        return [record["seq"] for record in candidate[len(candidate) - len(ops):]]

    def correct(self, target_seq: int, new_ops: Sequence[Mapping], note: str | None = None) -> list[int]:
        """Replace the step at ``target_seq`` and replay everything after it.

        Refused outright if any product version at or after ``target_seq`` has
        already been consumed by an analysis: old provenance is never silently
        rewritten.
        """
        if target_seq < 1 or target_seq > len(self.ops):
            raise LedgerError(f"cannot correct seq {target_seq}; log has {len(self.ops)} steps")
        # Snapshot in-memory state so any refusal below leaves it untouched.
        prev_ops = self.ops
        prev_state = (self.unknowns, self.containers, self.products)
        consumed = self._consumed_products()
        for p in self.products:
            if p.seq >= target_seq and p.product_id in consumed:
                raise LedgerError(
                    f"refusing correction: product {p.product_id} (container {p.container_id}, seq {p.seq}) "
                    f"has already been consumed by an analysis"
                )
        new_ops = list(new_ops)
        if not new_ops:
            raise LedgerError("correction must supply at least one replacement operation")
        prefix = self.ops[: target_seq - 1]
        suffix = self.ops[target_seq:]
        rebuilt: list[dict] = [dict(r) for r in prefix]
        seq = target_seq
        for op in new_ops:
            payload = json.loads(json.dumps(op))
            rebuilt.append({"seq": seq, "kind": op_kind(payload), "op": payload})
            seq += 1
        for record in suffix:
            rebuilt.append({"seq": seq, "kind": record["kind"], "op": json.loads(json.dumps(record["op"]))})
            seq += 1

        old_product_ids = {p.product_id for p in self.products}
        old_target_payload = json.dumps(self.ops[target_seq - 1]["op"], sort_keys=True)
        self._replay(rebuilt, persist=False)  # raises and leaves stored state untouched
        new_product_ids = {p.product_id for p in self.products}
        stale_consumed = (old_product_ids - new_product_ids) & self._consumed_products()
        if stale_consumed:
            self.unknowns, self.containers, self.products = prev_state
            raise LedgerError(f"refusing correction: consumed products would disappear: {sorted(stale_consumed)}")

        with self.conn:
            self.conn.execute("DELETE FROM ops")
            for record in rebuilt:
                self.conn.execute(
                    "INSERT INTO ops(seq, kind, payload) VALUES (?,?,?)",
                    (record["seq"], record["kind"], json.dumps(record["op"], sort_keys=True)),
                )
            self.conn.execute(
                "INSERT INTO corrections(target_seq, replaced_payload, new_payload, note, created_at) VALUES (?,?,?,?,?)",
                (
                    target_seq,
                    old_target_payload,
                    json.dumps([o for o in new_ops], sort_keys=True),
                    note,
                    datetime.now(timezone.utc).isoformat(timespec="seconds"),
                ),
            )
        self.ops = rebuilt
        self._replay(rebuilt, persist=False)
        self._persist_products(self.products)
        return [r["seq"] for r in rebuilt]

    def has_product(self, product_id: str) -> bool:
        return any(p.product_id == product_id for p in self.products)

    def mark_consumed(self, product_id: str, ref: str) -> None:
        """Record that ``product_id`` has been consumed by an analysis.

        Refuses unknown ids so a stale export cannot silently pin a product that
        no longer exists. Consumption makes the product version immutable: any
        later correction that would replace it is refused by :meth:`correct`.
        """
        self.product(product_id)  # raises if the id is not in the current ledger
        with self.conn:
            self.conn.execute(
                "INSERT OR IGNORE INTO consumed(product_id, ref, consumed_at) VALUES (?,?,?)",
                (product_id, ref, datetime.now(timezone.utc).isoformat(timespec="seconds")),
            )

    def consumed_products(self) -> list[dict]:
        rows = self.conn.execute("SELECT product_id, ref, consumed_at FROM consumed ORDER BY consumed_at, product_id").fetchall()
        return [dict(r) for r in rows]

    def _consumed_products(self) -> set[str]:
        try:
            rows = self.conn.execute("SELECT DISTINCT product_id FROM consumed").fetchall()
        except sqlite3.OperationalError:
            return set()
        return {r["product_id"] for r in rows}

    # ---------------------------------------------------------------- queries
    @property
    def latest_seq(self) -> int:
        return len(self.ops)

    def list_containers(self) -> list[Container]:
        return [self.containers[cid] for cid in sorted(self.containers)]

    def container(self, cid: str) -> Container:
        if cid not in self.containers:
            raise LedgerError(f"container {cid!r} does not exist")
        return self.containers[cid]

    def current_product(self, cid: str) -> Product:
        if cid not in self.containers:
            raise LedgerError(f"container {cid!r} does not exist")
        matches = [p for p in self.products if p.container_id == cid]
        if not matches:
            raise LedgerError(f"container {cid!r} has no product version")
        return max(matches, key=lambda p: p.version)

    def product(self, product_id: str) -> Product:
        for p in self.products:
            if p.product_id == product_id:
                return p
        raise LedgerError(f"product {product_id!r} not found in this ledger")

    def explain(self, cid: str) -> dict:
        c = self.container(cid)
        return {
            "container": c.cid,
            "kind": c.kind,
            "version": c.version,
            "volume": qstr(c.volume),
            "mass": c.mass.formatted(),
            "mass_terms": {k: qstr(v) for k, v in c.mass.sparse().items()},
            "concentration": (c.mass.scale(1 / c.volume)).formatted() if c.volume > 0 else "undefined (empty)",
            "origins": [
                {
                    "origin": origin,
                    "volume": qstr(vol),
                    "mass": mass.formatted(),
                    "mass_terms": {k: qstr(v) for k, v in mass.sparse().items()},
                }
                for origin, (vol, mass) in sorted(c.origins.items())
            ],
            "samples": list(self.unknowns),
        }

    def export_product(self, cid: str) -> dict:
        return self.current_product(cid).to_payload()

    def log(self) -> list[dict]:
        out = []
        for record in self.ops:
            out.append({"seq": record["seq"], "kind": record["kind"], "op": record["op"]})
        return out

    def corrections(self) -> list[dict]:
        rows = self.conn.execute("SELECT * FROM corrections ORDER BY id").fetchall()
        return [dict(r) for r in rows]


def product_from_payload(payload: Mapping) -> Product:
    """Parse an exported product payload.

    Only ``volume``, ``unknowns`` and ``mass`` are required; ``product_id`` is
    generated from the content when absent, so an external producer that does not
    know the ledger's id scheme can still submit a conforming product.
    """
    try:
        volume = q(payload["volume"], "volume")
        unknown_list = tuple(payload["unknowns"])
        mass_payload = payload["mass"]
        mass = LinExpr.from_rows(q(mass_payload["const"], "mass.const"), mass_payload.get("terms", {}), unknown_list)
        pid = payload.get("product_id")
        if not pid:
            digest = hashlib.sha256(
                json.dumps(
                    {"volume": qstr(volume), "unknowns": list(unknown_list),
                     "mass": {"const": qstr(mass.const), "terms": {k: qstr(v) for k, v in mass.sparse().items()}}},
                    sort_keys=True,
                ).encode()
            ).hexdigest()[:16]
            pid = f"ext-{digest}"
        origins = {}
        for item in payload.get("origins", []):
            mp = item["mass"]
            origins[item["origin"]] = (
                q(item["volume"], "origin.volume"),
                LinExpr.from_rows(q(mp["const"], "origin.mass.const"), mp.get("terms", {}), unknown_list),
            )
    except (KeyError, TypeError) as exc:
        raise LedgerError(f"malformed product payload: {exc}") from exc
    return Product(
        product_id=pid,
        container_id=payload.get("container", "?"),
        version=int(payload.get("version", 0)),
        seq=int(payload.get("seq", 0)),
        kind=payload.get("kind", "external"),
        volume=volume,
        mass=mass,
        origins=origins,
        unknowns=unknown_list,
    )
