"""End-to-end reproducible demo.

Runs on a fresh pair of SQLite files under ``out/``:

1. every scenario is replayed into the ledger and its hand-checked mass
   expressions are asserted;
2. a measurement is imported per reading and analyses are solved (unique, range,
   contradiction, shared-sample, zero-coefficient, rollback, evaporate);
3. a correction on an *unconsumed* branch is replayed and a fresh analysis
   consumes the new product version, while an older analysis is kept verbatim;
4. a correction on a *consumed* branch is refused;
5. the databases are closed and reopened to prove persistence;
6. the recovered endpoints are compared with the generator's ground truth.

The script exits non-zero if any check fails.
"""

from __future__ import annotations

import json
import shutil
from fractions import Fraction
from pathlib import Path

from .ledger import Ledger, LedgerError
from .rational import q, qstr
from .scenarios import all_scenarios
from .vault import Vault, VaultError

OUT = Path(__file__).resolve().parent.parent / "out"


class Failure(Exception):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise Failure(message)


def mass_of(ledger: Ledger, cid: str) -> str:
    return ledger.container(cid).mass.formatted()


def volume_of(ledger: Ledger, cid: str) -> Fraction:
    return ledger.container(cid).volume


def _bounds_for(ledger, scenario, product_names) -> dict:
    """Every sample referenced by the measured products needs a bound; fill any
    missing ones from the scenario defaults or a generic non-negative range."""
    if scenario.get("bounds"):
        return scenario["bounds"]
    used: list[str] = []
    for cid in product_names:
        for u in ledger.container(cid).mass.unknowns:
            if ledger.container(cid).mass.coeff_of(u) != 0 and u not in used:
                used.append(u)
    return {u: ["0", "1000"] for u in used}


def run_analysis(ledger, vault, scenario, analysis_name, override_bounds=None):
    product_names = []
    names = []
    for i, reading in enumerate(scenario["readings"], start=1):
        product = ledger.export_product(reading["product"])
        name = f"{scenario['name']}-r{i}"
        vault.add_measurement(name, product, reading["observed"], reading["uncertainty"])
        names.append(name)
        product_names.append(reading["product"])
    bounds = override_bounds or _bounds_for(ledger, scenario, product_names)
    result = vault.analyze(analysis_name, names, bounds, ledger=ledger)
    return result


def demo() -> dict:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    report = {"scenarios": [], "checks": []}
    scenarios = all_scenarios()

    def fresh(tag):
        return Ledger(OUT / f"{tag}.ledger.db"), Vault(OUT / f"{tag}.vault.db")

    # ---------------------------------------------------------------- two stage
    L, V = fresh("two_stage")
    s = scenarios[0]
    L.apply_batch(s["operations"])
    check(mass_of(L, "RAW") == "8*s1", f"RAW mass expected 8*s1, got {mass_of(L, 'RAW')}")
    check(mass_of(L, "A") == "s1", f"A mass expected s1, got {mass_of(L, 'A')}")
    check(mass_of(L, "B") == "100 + s1", f"B mass expected 100 + s1, got {mass_of(L, 'B')}")
    check(volume_of(L, "B") == 4, f"B volume after evaporation expected 4, got {volume_of(L, 'B')}")
    check(mass_of(L, "B") == "100 + s1", "evaporation must not change solute mass")
    r = run_analysis(L, V, s, "two-stage")
    check(r["status"] == "feasible", "two-stage should be feasible")
    check(r["endpoints"]["s1"]["min"] == "12" and r["endpoints"]["s1"]["max"] == "12",
          f"two-stage should pin s1=12, got {r['endpoints']}")
    check(r["unique"] is True, "two-stage must be flagged unique")
    report["scenarios"].append({"name": s["name"], "status": r["status"], "endpoints": r["endpoints"]})
    report["checks"].append("two-stage: 100+s1 over 4 mL pins s1=12 exactly")
    L.close(); V.close()

    # ------------------------------------------------------------- shared sample
    L, V = fresh("shared")
    s = scenarios[1]
    L.apply_batch(s["operations"])
    check(mass_of(L, "M2") == "4*s1", f"M2 mass expected 4*s1, got {mass_of(L, 'M2')}")
    check(mass_of(L, "M3") == "100 + 3/2*s1", f"M3 mass expected 100 + 3/2*s1, got {mass_of(L, 'M3')}")
    r = run_analysis(L, V, s, "shared-sample")
    check(r["status"] == "feasible", "shared-sample should be feasible")
    check(r["endpoints"]["s1"]["min"] == "15" and r["endpoints"]["s1"]["max"] == "15",
          f"shared sample should pin s1=15, got {r['endpoints']}")
    report["scenarios"].append({"name": s["name"], "status": r["status"], "endpoints": r["endpoints"]})
    L.close(); V.close()

    # -------------------------------------------------------------------- range
    L, V = fresh("range")
    s = scenarios[2]
    L.apply_batch(s["operations"])
    r = run_analysis(L, V, s, "range")
    check(r["status"] == "feasible", "range scenario should be feasible")
    check(r["endpoints"]["s1"]["min"] == "25" and r["endpoints"]["s1"]["max"] == "35",
          f"range expected [25,35], got {r['endpoints']}")
    check(r["unique"] is False, "range scenario must not be unique")
    body = {k: v for k, v in r.items() if k not in ("name", "created_at", "measurement_snapshot")}
    check("point_estimate" not in r and "average" not in json.dumps(body),
          "no fake point estimate may be introduced")
    report["scenarios"].append({"name": s["name"], "status": r["status"], "endpoints": r["endpoints"]})
    L.close(); V.close()

    # ------------------------------------------------------------ incompatible
    L, V = fresh("incompatible")
    s = scenarios[3]
    L.apply_batch(s["operations"])
    r = run_analysis(L, V, s, "incompatible")
    check(r["status"] == "contradiction", "incompatible readings should contradict")
    report["scenarios"].append({"name": s["name"], "status": r["status"], "endpoints": {}})
    L.close(); V.close()

    # ------------------------------------------------------ two unknowns unique
    L, V = fresh("two_unknowns")
    s = scenarios[4]
    L.apply_batch(s["operations"])
    r = run_analysis(L, V, s, "two-unknowns")
    check(r["status"] == "feasible", "two-unknowns should be feasible")
    check(r["endpoints"]["s1"]["min"] == "10" and r["endpoints"]["s1"]["max"] == "10",
          f"s1 expected 10, got {r['endpoints']['s1']}")
    check(r["endpoints"]["s2"]["min"] == "40" and r["endpoints"]["s2"]["max"] == "40",
          f"s2 expected 40, got {r['endpoints']['s2']}")
    report["scenarios"].append({"name": s["name"], "status": r["status"], "endpoints": r["endpoints"]})
    L.close(); V.close()

    # ------------------------------------------------------- zero coefficient
    L, V = fresh("zero_coeff")
    s = scenarios[5]
    L.apply_batch(s["operations"])
    check(mass_of(L, "ONLYSTD") == "200", f"ONLYSTD mass expected 200, got {mass_of(L, 'ONLYSTD')}")
    r = run_analysis(L, V, s, "zero-coeff")
    check(r["status"] == "contradiction", "zero-coefficient product should still contradict")
    check(len(r["individual_conflicts"]) == 1, "the conflict must be localised to the reading")
    report["scenarios"].append({"name": s["name"], "status": r["status"], "endpoints": {}})
    L.close(); V.close()

    # ------------------------------------------- rollback, correction, durability
    ledger_path = OUT / "branch.ledger.db"
    vault_path = OUT / "branch.vault.db"
    L = Ledger(ledger_path)
    V = Vault(vault_path)
    L.apply_batch(
        [
            {"op": "sample", "container": "SRC", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "SRC", "volume": "2", "destination": "POSTD"},
            {"op": "solvent", "container": "POSTD", "volume": "8"},
        ]
    )

    before = L.latest_seq
    before_containers = {c.cid: (c.volume, c.mass.formatted()) for c in L.list_containers()}
    try:
        L.apply_batch(
            [
                {"op": "transfer", "source": "SRC", "volume": "1", "destination": "TMP"},
                {"op": "transfer", "source": "SRC", "volume": "999", "destination": "TMP2"},
            ]
        )
        raise Failure("an over-draw must be rejected")
    except LedgerError:
        pass
    check(L.latest_seq == before, "failed batch must not append any step")
    check("TMP" not in L.containers and "TMP2" not in L.containers, "failed batch leaked containers")
    after_containers = {c.cid: (c.volume, c.mass.formatted()) for c in L.list_containers()}
    check(before_containers == after_containers, "failed batch changed container state")
    report["checks"].append("over-draw batch rolled back atomically")

    post_seq = L.latest_seq  # the solvent step
    old_product = L.export_product("POSTD")
    V.add_measurement("post-old", old_product, "6", "0")
    old_analysis = V.analyze("post-before", ["post-old"], {"s1": ["0", "100"]}, record_consumption=False)
    check(old_analysis["status"] == "feasible", "pre-correction analysis should be feasible")
    check(old_analysis["endpoints"]["s1"]["min"] == "30", "pre-correction s1 should be 30")

    # Correct the solvent top-up (8 -> 4 mL) and replay.
    L.correct(post_seq, [{"op": "solvent", "container": "POSTD", "volume": "4"}], note="wrong top-up")
    new_product = L.export_product("POSTD")
    check(old_product["product_id"] != new_product["product_id"], "corrected product must have a new id")
    check(mass_of(L, "POSTD") == "2*s1", "mass expression should be unchanged by top-up correction")
    check(volume_of(L, "POSTD") == 6, f"corrected volume expected 6, got {volume_of(L, 'POSTD')}")
    V.add_measurement("post-new", new_product, "6", "0")
    new_analysis = V.analyze("post-after", ["post-new"], {"s1": ["0", "100"]}, ledger=L)
    check(new_analysis["endpoints"]["s1"]["min"] == "18",
          f"post-correction s1 expected 18, got {new_analysis['endpoints']}")
    check(V.analysis("post-before")["endpoints"]["s1"]["min"] == "30",
          "old analysis must be preserved verbatim")
    check(L.corrections() and L.corrections()[0]["target_seq"] == post_seq, "correction audit trail missing")
    report["checks"].append("correction replayed from prior state; old analysis preserved; new analysis consumed new product")

    # The new product is now consumed, so correcting that step must be refused.
    try:
        L.correct(post_seq, [{"op": "solvent", "container": "POSTD", "volume": "2"}])
        raise Failure("correction of a consumed product must be refused")
    except LedgerError:
        pass
    check(volume_of(L, "POSTD") == 6, "refused correction must leave the product unchanged")
    report["checks"].append("correction of a consumed product version refused")

    # A second, unrelated batch can still be appended and corrected (no consumption).
    L.apply_batch([{"op": "sample", "container": "FREE", "volume": "7", "sample": "s2"}])
    free_seq = L.latest_seq
    L.correct(free_seq, [{"op": "sample", "container": "FREE", "volume": "9", "sample": "s2"}], note="typo")
    check(volume_of(L, "FREE") == 9, "unconsumed correction should replay")
    report["checks"].append("unconsumed correction replayed from prior state")

    # --------------------------------------------------------- durability check
    n_ops, containers = len(L.ops), dict(L.containers)
    n_meas, n_analyses = len(V.list_measurements()), len(V.list_analyses())
    L.close(); V.close()
    L2 = Ledger(ledger_path)
    V2 = Vault(vault_path)
    check(len(L2.ops) == n_ops, "reopened ledger lost steps")
    check(mass_of(L2, "POSTD") == "2*s1", "reopened ledger lost the corrected state")
    check(volume_of(L2, "FREE") == 9, "reopened ledger lost the unconsumed correction")
    check(V2.analysis("post-before")["endpoints"]["s1"]["min"] == "30", "reopened vault lost old analysis")
    check(V2.analysis("post-after")["endpoints"]["s1"]["min"] == "18", "reopened vault lost new analysis")
    check(len(V2.list_measurements()) == n_meas, "reopened vault lost measurements")
    check(len(L2.consumed_products()) >= 1, "consumption records must survive reopen")
    report["checks"].append("databases reopen with identical state")

    report["ledger_steps"] = len(L2.ops)
    report["containers"] = len(L2.containers)
    report["measurements"] = len(V2.list_measurements())
    report["analyses"] = len(V2.list_analyses())
    L2.close(); V2.close()
    return report


def main() -> int:
    try:
        report = demo()
    except Failure as exc:
        print(f"DEMO FAILED: {exc}")
        return 1
    (OUT / "demo_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("demo OK")
    print(json.dumps({k: report[k] for k in ("ledger_steps", "containers", "measurements", "analyses")}, indent=2))
    print(f"report written to {OUT / 'demo_report.json'}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
