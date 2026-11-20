"""Synthetic scenario generator.

Every scenario is a plain JSON operation batch (plus a vault plan). The generator
also reports the *true* sample concentrations it used, so the demo can compare
the exact recovered intervals with ground truth — the solver never sees them.

Scenarios are deliberately hand-checkable (small integer volumes and masses).
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from .rational import q, qstr


def _op(**payload):
    return payload


def scenario_two_stage() -> dict:
    """Two-stage dilution of one sample with a standard spike and evaporation.

    Hand check with s1 = 12 mg/mL (unknown to the solver):
      RAW 10 mL -> 120 mg.
      A = 2 mL of RAW + 8 mL solvent -> 24 mg in 10 mL = 2.4 mg/mL.
      B = 5 mL of A + 1 mL of STD(100) -> 12 mg + 100 mg = 112 mg in 6 mL.
      evaporate B to 4 mL -> still 112 mg in 4 mL = 28 mg/mL.
      The flask therefore holds ``100 + s1`` mg in 4 mL, i.e. concentration
      ``25 + s1/4``. A reading of 28 pins s1 = 12 exactly.

    The recovered solution also demonstrates that multiplying the reading by the
    last dilution factor (4/6) would be wrong: 28 * 6/4 = 42 is not the sample
    concentration because 100 mg of standard is in the same flask.
    """
    ops = [
        _op(op="standard", container="STD", volume="10", concentration="100", label="100 mg/mL standard"),
        _op(op="sample", container="RAW", volume="10", sample="s1", label="raw sample"),
        _op(op="transfer", source="RAW", volume="2", destination="A"),
        _op(op="solvent", container="A", volume="8"),
        _op(op="mix", sources=[{"container": "A", "volume": "5"}, {"container": "STD", "volume": "1"}], destination="B"),
        _op(op="evaporate", container="B", to="4"),
    ]
    truth = {"s1": "12"}
    readings = [{"product": "B", "observed": "28", "uncertainty": "0", "bounds": {"s1": ["0", "100"]}}]
    return {"name": "two_stage_spike_evap", "operations": ops, "truth": truth, "readings": readings}


def scenario_shared_sample() -> dict:
    """One sample feeding two *different* final mixtures.

    s1 = 15 mg/mL.
      RAW 8 mL -> 120 mg.
      M1 = 2 mL RAW + 2 mL solvent -> 30 mg in 4 mL  (7.5 mg/mL).
      M2 = 4 mL RAW + 4 mL solvent -> 60 mg in 8 mL  (7.5 mg/mL).
    Same concentration twice but *independently* derived; both must be included.
    A third mix contaminates a standard:
      M3 = 3 mL M1 + 2 mL STD(50) -> 22.5 + 100 = 122.5 mg in 5 mL = 24.5 mg/mL.
    Readings of 7.5 and 24.5 over-determine s1; if one is perturbed the joint
    solution becomes a narrow range or a contradiction.
    """
    ops = [
        _op(op="standard", container="STD", volume="10", concentration="50"),
        _op(op="sample", container="RAW", volume="8", sample="s1"),
        _op(op="transfer", source="RAW", volume="2", destination="M1"),
        _op(op="solvent", container="M1", volume="2"),
        _op(op="transfer", source="RAW", volume="4", destination="M2"),
        _op(op="solvent", container="M2", volume="4"),
        _op(op="mix", sources=[{"container": "M1", "volume": "3"}, {"container": "STD", "volume": "2"}], destination="M3"),
    ]
    truth = {"s1": "15"}
    readings = [
        # M1 has 30 mg in 4 mL -> 7.5; M2 has 60 mg in 8 mL -> 7.5
        {"product": "M1", "observed": "15/2", "uncertainty": "0"},
        {"product": "M3", "observed": "49/2", "uncertainty": "0"},
    ]
    return {
        "name": "shared_sample_two_mixtures",
        "operations": ops,
        "truth": truth,
        "readings": readings,
        "bounds": {"s1": ["0", "100"]},
    }


def scenario_range() -> dict:
    """A single interval reading leaves a genuine range, never a fake average.

    s1 is not pinned: D = 2 mL RAW + 8 mL solvent, 10 mL.
    Concentration = s1/5. Reading 6 ± 1 => s1 ∈ [25, 35].
    """
    ops = [
        _op(op="sample", container="RAW", volume="10", sample="s1"),
        _op(op="transfer", source="RAW", volume="2", destination="D"),
        _op(op="solvent", container="D", volume="8"),
    ]
    return {
        "name": "range_from_uncertainty",
        "operations": ops,
        "truth": {"s1": "30"},
        "readings": [{"product": "D", "observed": "6", "uncertainty": "1"}],
        "bounds": {"s1": ["0", "100"]},
    }


def scenario_incompatible() -> dict:
    """Two readings of the same product that cannot both be true."""
    ops = [
        _op(op="sample", container="RAW", volume="10", sample="s1"),
        _op(op="transfer", source="RAW", volume="2", destination="D"),
        _op(op="solvent", container="D", volume="8"),
    ]
    return {
        "name": "incompatible_readings",
        "operations": ops,
        "truth": {"s1": "30"},
        "readings": [
            {"product": "D", "observed": "6", "uncertainty": "0"},
            {"product": "D", "observed": "7", "uncertainty": "0"},
        ],
        "bounds": {"s1": ["0", "100"]},
    }


def scenario_two_unknowns_unique() -> dict:
    """Two samples, two independent linear equations, exact unique solution.

    s1 = 10 (RAW1), s2 = 40 (RAW2).
      F1 = 1 mL RAW1 + 1 mL solvent  -> mass s1, in 2 mL.
      F2 = 2 mL RAW2 + 2 mL solvent  -> mass 2*s2, in 4 mL.
      G  = 1 mL F1 + 1 mL F2 -> (s1 + 2*s2)/... careful: F1 has mass s1 mg in
           2 mL, so 1 mL contributes s1/2; F2 1 mL contributes (2 s2)/4 = s2/2.
           G mass = (s1+s2)/2, volume 2.
    Instead use pure transfers to keep coefficients integral:
      P1 = 2 mL RAW1 diluted to 4 mL -> concentration s1/2.
      P2 = 3 mL RAW2 diluted to 6 mL -> concentration s2/2.
    Readings 5 and 20 pin s1 = 10, s2 = 40.
    """
    ops = [
        _op(op="sample", container="R1", volume="10", sample="s1"),
        _op(op="sample", container="R2", volume="10", sample="s2"),
        _op(op="transfer", source="R1", volume="2", destination="P1"),
        _op(op="solvent", container="P1", volume="2"),
        _op(op="transfer", source="R2", volume="3", destination="P2"),
        _op(op="solvent", container="P2", volume="3"),
    ]
    return {
        "name": "two_unknowns_unique",
        "operations": ops,
        "truth": {"s1": "10", "s2": "40"},
        "readings": [
            {"product": "P1", "observed": "5", "uncertainty": "0"},
            {"product": "P2", "observed": "20", "uncertainty": "0"},
        ],
        "bounds": {"s1": ["0", "100"], "s2": ["0", "100"]},
    }


def scenario_zero_coefficient() -> dict:
    """A product with *no* sample term can still contradict a reading.

    STD alone (100 mg/mL), but the reading says 50 => contradiction that no
    sample choice can fix. Demonstrates that zero-coefficient rows are real.
    """
    ops = [
        _op(op="standard", container="STD", volume="5", concentration="100"),
        _op(op="sample", container="RAW", volume="5", sample="s1"),
        _op(op="mix", sources=[{"container": "STD", "volume": "2"}], destination="ONLYSTD"),
    ]
    return {
        "name": "zero_coefficient_contradiction",
        "operations": ops,
        "truth": {"s1": "10"},
        "readings": [{"product": "ONLYSTD", "observed": "50", "uncertainty": "0"}],
        "bounds": {"s1": ["0", "100"]},
    }


def all_scenarios() -> list[dict]:
    return [
        scenario_two_stage(),
        scenario_shared_sample(),
        scenario_range(),
        scenario_incompatible(),
        scenario_two_unknowns_unique(),
        scenario_zero_coefficient(),
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="prepquant.scenarios")
    parser.add_argument("--out", default="out/scenarios.json")
    args = parser.parse_args(argv)
    data = all_scenarios()
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {len(data)} scenarios to {target}")
    for s in data:
        print(f"  {s['name']}: {len(s['operations'])} operations, {len(s['readings'])} reading(s)")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
