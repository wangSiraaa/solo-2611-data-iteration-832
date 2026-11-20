import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_cli(tmp_path, *args, input_text=None):
    cmd = [
        sys.executable,
        "-m",
        "prepquant.cli",
        "--ledger",
        str(tmp_path / "l.db"),
        "--vault",
        str(tmp_path / "v.db"),
        *args,
    ]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, input=input_text)


def test_cli_full_flow(tmp_path):
    assert run_cli(tmp_path, "init").returncode == 0

    ops = [
        {"op": "standard", "container": "STD", "volume": "10", "concentration": "100"},
        {"op": "sample", "container": "RAW", "volume": "10", "sample": "s1"},
        {"op": "transfer", "source": "RAW", "volume": "2", "destination": "D"},
        {"op": "solvent", "container": "D", "volume": "8"},
        {"op": "mix", "sources": [{"container": "D", "volume": "5"}, {"container": "STD", "volume": "1"}], "destination": "M"},
    ]
    ops_file = tmp_path / "ops.json"
    ops_file.write_text(json.dumps(ops), encoding="utf-8")
    r = run_cli(tmp_path, "prep", str(ops_file))
    assert r.returncode == 0, r.stderr

    r = run_cli(tmp_path, "explain", "M", "--json")
    assert r.returncode == 0, r.stderr
    explain = json.loads(r.stdout)
    assert explain["volume"] == "6"
    assert explain["mass"] == "100 + s1"
    assert {o["origin"] for o in explain["origins"]} == {"STD", "RAW", "solvent"}

    out = tmp_path / "m.json"
    r = run_cli(tmp_path, "export", "M", "--out", str(out))
    assert r.returncode == 0, r.stderr
    payload = json.loads(out.read_text(encoding="utf-8"))
    # reading 25 mg/mL at 6 mL -> 150 mg = 100 + s1 -> s1 = 50
    r = run_cli(tmp_path, "measure", "--name", "r1", "--product-file", str(out),
                "--observed", "25", "--uncertainty", "0")
    assert r.returncode == 0, r.stderr
    r = run_cli(tmp_path, "analyze", "a1", "--measurement", "r1", "--bound", "s1:0:100")
    assert r.returncode == 0, r.stderr
    assert "50" in r.stdout
    assert "unique" in r.stdout

    r = run_cli(tmp_path, "show-analysis", "a1")
    assert r.returncode == 0, r.stderr
    assert "status: feasible" in r.stdout

    r = run_cli(tmp_path, "consumed")
    assert payload["product_id"] in r.stdout

    # reopening neither loses the step log nor the analysis
    r = run_cli(tmp_path, "init")
    assert "5 steps" in r.stdout
    assert run_cli(tmp_path, "show-analysis", "a1").returncode == 0


def test_cli_rejects_bad_batch_without_partial_apply(tmp_path):
    ops = [
        {"op": "sample", "container": "R", "volume": "2", "sample": "s1"},
        {"op": "transfer", "source": "R", "volume": "5", "destination": "X"},
    ]
    ops_file = tmp_path / "bad.json"
    ops_file.write_text(json.dumps(ops), encoding="utf-8")
    assert run_cli(tmp_path, "prep", str(ops_file)).returncode != 0
    r = run_cli(tmp_path, "containers")
    assert json.loads(r.stdout) == []
    r = run_cli(tmp_path, "log")
    assert r.stdout.strip() == ""


def test_cli_correct(tmp_path):
    ops = [
        {"op": "sample", "container": "R", "volume": "10", "sample": "s1"},
        {"op": "transfer", "source": "R", "volume": "2", "destination": "D"},
        {"op": "solvent", "container": "D", "volume": "8"},
    ]
    f = tmp_path / "ops.json"
    f.write_text(json.dumps(ops), encoding="utf-8")
    assert run_cli(tmp_path, "prep", str(f)).returncode == 0
    repl = tmp_path / "fix.json"
    repl.write_text(json.dumps([{"op": "solvent", "container": "D", "volume": "4"}]), encoding="utf-8")
    r = run_cli(tmp_path, "correct", "3", str(repl), "--note", "typo")
    assert r.returncode == 0, r.stderr
    r = run_cli(tmp_path, "explain", "D", "--json")
    assert json.loads(r.stdout)["volume"] == "6"
    corrections = json.loads(run_cli(tmp_path, "corrections").stdout)
    assert corrections[0]["target_seq"] == 3
