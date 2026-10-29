"""qpAdm's CLI output is complete (TASKS.md H4).

Weights with SEs and z in text, TSV and JSON; the rank tests (rankdrop)
and popdrop; the blocks and SNPs used; p-values and SEs at full precision
in TSV and JSON and never as 0.0000 in text; and rotations as one TSV row
per model. The numbers are those of the Python API, checked against
admixtools2 in tests/golden.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pandas as pd
import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
DG6 = ROOT / "tests" / "golden" / "admixtools2" / "v2.0.10" / "dg_6pop" / "f2_dir"
MODEL = ["--target", "French.DG", "--left", "Han.DG,Yoruba.DG", "--right", "Mbuti.DG,Chimp.REF,Papuan.DG"]


@pytest.fixture(scope="module")
def f2_dir(tmp_path_factory) -> str:
    out = tmp_path_factory.mktemp("dg6") / "f2"
    pm.read_admixtools2_f2(DG6).save(str(out))
    return str(out)


@pytest.fixture(scope="module")
def result(f2_dir) -> pm.QpAdmResult:
    return pm.qpadm(pm.read_f2(f2_dir), ["Han.DG", "Yoruba.DG"], ["Mbuti.DG", "Chimp.REF", "Papuan.DG"],
                    "French.DG")


def popmix(*args: str) -> subprocess.CompletedProcess:
    r = subprocess.run(["popmix", *args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


def test_the_result_has_z_rankdrop_and_counts(result, f2_dir):
    f2 = pm.read_f2(f2_dir)
    assert result.weight_z == tuple(w / se for w, se in zip(result.weights, result.weight_se))
    assert [t.rank for t in result.rankdrop] == [0, 1]
    assert result.n_blocks == f2.n_blocks
    assert result.n_snps == sum(f2.block_lengths)


def test_text_shows_z_counts_rank_tests_and_small_p(result, f2_dir):
    out = popmix("qpadm", "--f2-dir", f2_dir, *MODEL).stdout
    assert f"{result.n_blocks} jackknife blocks, {result.n_snps} SNPs" in out
    han = next(line for line in out.splitlines() if line.strip().startswith("Han.DG"))
    assert han.split() == ["Han.DG", f"{result.weights[0]:.4f}", f"{result.weight_se[0]:.3g}",
                           f"{result.weight_z[0]:.1f}"]
    assert "Rank tests (rankdrop):" in out
    rank0 = result.rankdrop[0]
    assert f"{rank0.chisq:.4f}" in out
    assert "0.0000" not in out


def test_tsv_has_full_precision_and_z(result, f2_dir):
    out = popmix("qpadm", "--f2-dir", f2_dir, *MODEL, "--format", "tsv").stdout
    df = pd.read_csv(pd.io.common.StringIO(out), sep="\t")
    assert list(df.columns) == ["target", "source", "weight", "se", "z", "chisq", "dof", "p", "feasible",
                                "n_blocks", "n_snps"]
    assert df.weight.tolist() == pytest.approx(result.weights, rel=1e-11)
    assert df.se.tolist() == pytest.approx(result.weight_se, rel=1e-11)
    assert df.p[0] == pytest.approx(result.p, rel=1e-11)


def test_json_has_everything(result, f2_dir):
    d = json.loads(popmix("qpadm", "--f2-dir", f2_dir, *MODEL, "--format", "json").stdout)
    assert d["z"] == dict(zip(result.sources, result.weight_z))
    assert d["rankdrop"] == [{"f4rank": t.rank, "dof": t.dof, "chisq": t.chisq, "p": t.p} for t in result.rankdrop]
    assert (d["n_blocks"], d["n_snps"]) == (result.n_blocks, result.n_snps)


def test_an_output_prefix_writes_every_table(result, f2_dir, tmp_path):
    popmix("qpadm", "--f2-dir", f2_dir, *MODEL, "-o", str(tmp_path / "m"))
    weights = pd.read_csv(tmp_path / "m.weights.tsv", sep="\t")
    rankdrop = pd.read_csv(tmp_path / "m.rankdrop.tsv", sep="\t")
    popdrop = pd.read_csv(tmp_path / "m.popdrop.tsv", sep="\t", dtype={"pattern": str})
    assert weights.source.tolist() == list(result.sources)
    assert rankdrop.f4rank.tolist() == [0, 1]
    assert popdrop.pattern.tolist() == [r.pattern for r in result.popdrop]
    assert popdrop.p.tolist() == pytest.approx([r.p for r in result.popdrop], rel=1e-11)
    assert json.loads((tmp_path / "m.json").read_text(encoding="utf-8"))["p"] == pytest.approx(result.p, rel=1e-15)


def test_rotation_writes_one_row_per_model(f2_dir):
    out = popmix("qpadm-rotate", "--f2-dir", f2_dir, "--target", "French.DG", "--leftright",
                 "Han.DG,Yoruba.DG,Papuan.DG", "--rightfix", "Mbuti.DG,Chimp.REF").stdout
    df = pd.read_csv(pd.io.common.StringIO(out), sep="\t")
    f2 = pm.read_f2(f2_dir)
    want = pm.qpadm_rotate(f2, ["Han.DG", "Yoruba.DG", "Papuan.DG"], "French.DG", ["Mbuti.DG", "Chimp.REF"],
                           full_results=True)
    assert len(df) == len(want)
    assert df.left.tolist() == [",".join(r.sources) for r in want]
    assert df.p.tolist() == pytest.approx([r.p for r in want], rel=1e-11)
    first = want[0]
    assert df.weights[0] == ",".join(f"{w:.12g}" for w in first.weights)
    assert df.status[0] == first.assess().status
    assert math.isfinite(df.chisq[0])
