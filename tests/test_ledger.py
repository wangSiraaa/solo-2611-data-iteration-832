import os
import tempfile
from fractions import Fraction

import pytest

from prepquant.ledger import Ledger, LedgerError, product_from_payload


@pytest.fixture()
def ledger():
    d = tempfile.mkdtemp()
    L = Ledger(os.path.join(d, "ledger.db"))
    yield L
    L.close()


def op(**kw):
    return kw


def test_transfer_splits_pro_ratio(ledger):
    ledger.apply_batch(
        [
            op(op="standard", container="S", volume="10", concentration="80"),
            op(op="transfer", source="S", volume="3", destination="D"),
        ]
    )
    assert ledger.container("S").volume == 7
    assert ledger.container("S").mass.const == 560
    assert ledger.container("D").volume == 3
    assert ledger.container("D").mass.const == 240
    # provenance books total the container exactly
    for cid in ("S", "D"):
        c = ledger.container(cid)
        assert sum(v for v, _ in c.origins.values()) == c.volume
        total_mass = sum((m.const for _, m in c.origins.values()), Fraction(0))
        assert total_mass == c.mass.const


def test_cannot_use_same_liquid_twice(ledger):
    ledger.apply_batch([op(op="sample", container="R", volume="5", sample="s1")])
    with pytest.raises(LedgerError):
        ledger.apply_batch(
            [
                op(op="transfer", source="R", volume="3", destination="A"),
                op(op="transfer", source="R", volume="3", destination="B"),
            ]
        )
    # whole batch rolled back
    assert ledger.latest_seq == 1
    assert "A" not in ledger.containers and "B" not in ledger.containers
    assert ledger.container("R").volume == 5


def test_mix_and_evaporate(ledger):
    ledger.apply_batch(
        [
            op(op="sample", container="R", volume="10", sample="s1"),
            op(op="mix", sources=[{"container": "R", "volume": "4"}], destination="M"),
            op(op="solvent", container="M", volume="4"),
            op(op="evaporate", container="M", to="2"),
        ]
    )
    c = ledger.container("M")
    assert c.volume == 2
    assert c.mass.coeffs[c.mass.unknowns.index("s1")] == 4
    explain = ledger.explain("M")
    assert explain["mass"] == "4*s1"


def test_evaporation_must_leave_positive_volume(ledger):
    ledger.apply_batch([op(op="standard", container="S", volume="5", concentration="10")])
    with pytest.raises(LedgerError):
        ledger.apply_batch([op(op="evaporate", container="S", remove="5")])
    with pytest.raises(LedgerError):
        ledger.apply_batch([op(op="evaporate", container="S", to="0")])
    with pytest.raises(LedgerError):
        ledger.apply_batch([op(op="evaporate", container="S", to="6")])
    assert ledger.container("S").volume == 5


def test_discard_then_reuse_identity(ledger):
    ledger.apply_batch(
        [
            op(op="standard", container="S", volume="5", concentration="10"),
            op(op="discard", container="S"),
            op(op="standard", container="S", volume="2", concentration="20"),
        ]
    )
    assert ledger.container("S").volume == 2
    assert ledger.container("S").mass.const == 40


def test_origin_book_reports_removed_material(ledger):
    ledger.apply_batch(
        [
            op(op="standard", container="S", volume="10", concentration="10"),
            op(op="mix", sources=[{"container": "S", "volume": "4"}], destination="M"),
            op(op="evaporate", container="M", to="2"),
        ]
    )
    e = ledger.explain("M")
    origins = {o["origin"]: (o["volume"], o["mass"]) for o in e["origins"]}
    # 4 mL drawn, 40 mg solute; evaporation halves volume but keeps all 40 mg
    assert origins["S"] == ("2", "40")
    assert Fraction(e["volume"]) == 2


def test_product_versions_are_content_addressed(ledger):
    ledger.apply_batch([op(op="standard", container="S", volume="10", concentration="10")])
    first = ledger.export_product("S")
    ledger.apply_batch([op(op="solvent", container="S", volume="10")])
    second = ledger.export_product("S")
    assert first["product_id"] != second["product_id"]
    assert first["mass"] == second["mass"]  # dilution adds no solute
    assert ledger.product(first["product_id"]).volume == 10


def test_correction_replays_later_steps(ledger):
    ledger.apply_batch(
        [
            op(op="sample", container="R", volume="10", sample="s1"),
            op(op="transfer", source="R", volume="2", destination="D"),
            op(op="solvent", container="D", volume="8"),
        ]
    )
    ledger.correct(2, [op(op="transfer", source="R", volume="3", destination="D")])
    # later solvent step was replayed, under the same seq
    assert ledger.latest_seq == 3
    assert ledger.container("D").volume == 11
    assert ledger.container("D").mass.formatted() == "3*s1"
    assert ledger.corrections()[0]["target_seq"] == 2


def test_correction_failure_leaves_state(ledger):
    ledger.apply_batch(
        [
            op(op="sample", container="R", volume="10", sample="s1"),
            op(op="transfer", source="R", volume="2", destination="D"),
        ]
    )
    snapshot = {c.cid: (c.volume, c.mass.formatted()) for c in ledger.list_containers()}
    with pytest.raises(LedgerError):
        ledger.correct(1, [op(op="sample", container="X", volume="1", sample="s2"), op(op="transfer", source="X", volume="99", destination="Y")])
    assert ledger.latest_seq == 2
    assert {c.cid: (c.volume, c.mass.formatted()) for c in ledger.list_containers()} == snapshot


def test_consumed_product_blocks_correction(ledger):
    ledger.apply_batch(
        [
            op(op="sample", container="R", volume="10", sample="s1"),
            op(op="transfer", source="R", volume="2", destination="D"),
        ]
    )
    pid = ledger.export_product("D")["product_id"]
    ledger.mark_consumed(pid, ref="test")
    with pytest.raises(LedgerError):
        ledger.correct(2, [op(op="transfer", source="R", volume="3", destination="D")])
    assert ledger.container("D").volume == 2
    assert ledger.consumed_products()[0]["product_id"] == pid


def test_unknown_container_and_bad_fields(ledger):
    with pytest.raises(LedgerError):
        ledger.apply_batch([op(op="transfer", source="nope", volume="1", destination="D")])
    with pytest.raises(LedgerError):
        ledger.apply_batch([op(op="standard", container="S", volume="1", concentration="1", surprise=1)])
    with pytest.raises(LedgerError):
        ledger.apply_batch([op(op="standard", container="S", volume="-1", concentration="1")])
    with pytest.raises(LedgerError):
        ledger.apply_batch([op(op="standard", container="S", volume="0", concentration="1")])


def test_sample_limit(ledger):
    ops = [op(op="sample", container=f"R{i}", volume="1", sample=f"s{i}") for i in range(7)]
    with pytest.raises(LedgerError):
        ledger.apply_batch(ops)


def test_external_product_payload_without_id():
    prod = product_from_payload({"volume": "4", "unknowns": ["a"], "mass": {"const": "2", "terms": {"a": "3"}}})
    assert prod.product_id.startswith("ext-")
    assert prod.mass.formatted() == "2 + 3*a"
