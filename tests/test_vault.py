import os
import tempfile

import pytest

from prepquant.ledger import Ledger
from prepquant.vault import Vault, VaultError


@pytest.fixture()
def env():
    d = tempfile.mkdtemp()
    L = Ledger(os.path.join(d, "ledger.db"))
    V = Vault(os.path.join(d, "vault.db"))
    yield L, V
    L.close()
    V.close()


def setup_product(L, ops):
    L.apply_batch(ops)
    return L.export_product("M")


def test_unique_exact_solution(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "8"},
        ]
    )
    V.add_measurement("r", L.export_product("M"), "6", "0")
    res = V.analyze("a", ["r"], {"s1": ["0", "100"]})
    assert res["status"] == "feasible"
    assert res["unique"] is True
    assert res["endpoints"]["s1"]["min"] == "30"
    assert res["endpoints"]["s1"]["max"] == "30"


def test_range_is_not_averaged(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "8"},
        ]
    )
    V.add_measurement("r", L.export_product("M"), "6", "1")
    res = V.analyze("a", ["r"], {"s1": ["0", "100"]})
    assert res["endpoints"]["s1"] == {"min": "25", "max": "35", "unique": False}
    assert res["unique"] is False


def test_uncertainty_clipped_at_zero(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "8"},
        ]
    )
    # interval [-2, 2] intersected with [0, inf) -> [0,2]; concentration = s1/5
    V.add_measurement("r", L.export_product("M"), "0", "2")
    res = V.analyze("a", ["r"], {"s1": ["0", "100"]})
    assert res["endpoints"]["s1"] == {"min": "0", "max": "10", "unique": False}


def test_sample_bounds_are_applied(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "8"},
        ]
    )
    V.add_measurement("r", L.export_product("M"), "6", "0")  # -> s1 = 30
    res = V.analyze("a", ["r"], {"s1": ["0", "20"]})
    assert res["status"] == "contradiction"


def test_shared_sample_joint_solve(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "sample", "container": "Q", "volume": "10", "sample": "s2"},
            {"op": "mix", "sources": [{"container": "R", "volume": "2"}, {"container": "Q", "volume": "2"}], "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "6"},
        ]
    )
    V.add_measurement("r", L.export_product("M"), "6", "0")
    res = V.analyze("a", ["r"], {"s1": ["0", "100"], "s2": ["0", "100"]})
    assert res["witness_point"]
    assert res["endpoints"]["s1"]["max"] == "30"
    assert res["endpoints"]["s2"]["max"] == "30"
    assert res["unique"] is False


def test_multiple_measurements_of_shared_sample(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "sample", "container": "Q", "volume": "10", "sample": "s2"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "P"},
            {"op": "solvent", "container": "P", "volume": "8"},
            {"op": "mix", "sources": [{"container": "R", "volume": "2"}, {"container": "Q", "volume": "2"}], "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "6"},
        ]
    )
    V.add_measurement("p", L.export_product("P"), "6", "0")  # s1 = 30
    V.add_measurement("m", L.export_product("M"), "14", "0")  # (2*s1+2*s2)/10 = 14 -> s2=40
    res = V.analyze("a", ["p", "m"], {"s1": ["0", "100"], "s2": ["0", "100"]})
    assert res["endpoints"]["s1"]["min"] == "30" and res["endpoints"]["s2"]["min"] == "40"
    assert res["unique"] is True


def test_contradiction_and_zero_coefficient(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "standard", "container": "S", "volume": "5", "concentration": "100"},
            {"op": "sample", "container": "R", "volume": "5", "sample": "s1"},
            {"op": "mix", "sources": [{"container": "S", "volume": "2"}], "destination": "M"},
        ]
    )
    V.add_measurement("r", L.export_product("M"), "50", "0")
    res = V.analyze("a", ["r"], {"s1": ["0", "100"]})
    assert res["status"] == "contradiction"
    assert res["individual_conflicts"][0]["measurement"] == "r"


def test_external_product_format(env):
    L, V = env
    external = {
        "volume": "4",
        "unknowns": ["s1"],
        "mass": {"const": "8", "terms": {"s1": "2"}},
    }
    V.add_measurement("ext", external, "10", "0")  # (8+2 s1)/4 = 10 -> s1 = 16
    res = V.analyze("a", ["ext"], {"s1": ["0", "100"]})
    assert res["endpoints"]["s1"]["min"] == "16"
    assert res["endpoints"]["s1"]["max"] == "16"


def test_measurement_immutability_and_missing_bounds(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "M"},
        ]
    )
    V.add_measurement("r", L.export_product("M"), "2", "0")
    with pytest.raises(VaultError):
        V.add_measurement("r", L.export_product("M"), "3", "0")
    with pytest.raises(VaultError):
        V.analyze("a", ["r"], {})
    with pytest.raises(VaultError):
        V.analyze("a", ["r"], {"s1": ["10", "5"]})
    with pytest.raises(VaultError):
        V.analyze("a", ["missing"], {"s1": ["0", "5"]})


def test_reanalysis_keeps_old_result_and_consumes_new(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "8"},
        ]
    )
    p1 = L.export_product("M")
    V.add_measurement("r1", p1, "6", "0")
    a1 = V.analyze("first", ["r1"], {"s1": ["0", "100"]}, ledger=L)
    assert a1["endpoints"]["s1"]["min"] == "30"
    assert L.has_product(p1["product_id"]) and L.consumed_products()

    # Append an unrelated step; the old analysis must not change.
    L.apply_batch([{"op": "standard", "container": "S", "volume": "1", "concentration": "1"}])
    assert V.analysis("first")["endpoints"]["s1"]["min"] == "30"

    # A second analysis against the same product is allowed and also consumes it
    # (repeat consumption is idempotent per ref).
    a2 = V.analyze("second", ["r1"], {"s1": ["0", "100"]}, ledger=L)
    assert a2["endpoints"]["s1"]["min"] == "30"
    assert V.analysis("first")["endpoints"] == a2["endpoints"]


def test_analysis_snapshot_survives_correction_of_unconsumed_branch(env):
    L, V = env
    L.apply_batch(
        [
            {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
            {"op": "transfer", "source": "R", "volume": "2", "destination": "M"},
            {"op": "solvent", "container": "M", "volume": "8"},
        ]
    )
    old = L.export_product("M")
    V.add_measurement("old", old, "6", "0")
    a = V.analyze("old-analysis", ["old"], {"s1": ["0", "100"]}, record_consumption=False)
    assert a["endpoints"]["s1"]["min"] == "30"

    L.correct(3, [{"op": "solvent", "container": "M", "volume": "4"}])
    new = L.export_product("M")
    assert new["product_id"] != old["product_id"]
    V.add_measurement("new", new, "6", "0")
    b = V.analyze("new-analysis", ["new"], {"s1": ["0", "100"]}, ledger=L)
    assert b["endpoints"]["s1"]["min"] == "18"
    # Old analysis is stored verbatim, with its own expression snapshot.
    assert V.analysis("old-analysis")["endpoints"]["s1"]["min"] == "30"
    assert V.analysis("old-analysis")["measurement_snapshot"][0]["payload"]["product_id"] == old["product_id"]
