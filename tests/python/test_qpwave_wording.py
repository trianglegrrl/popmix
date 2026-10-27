"""qpWave says what its rank shows, at a level the user chooses (TASKS.md G2).

The estimated rank is the lowest rank not rejected at ``alpha``: a lower
bound on the number of streams of ancestry minus one, limited by the right
populations and the data. Not rejecting rank 0 does not show a clade.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
DG6 = str(ROOT / "tests" / "golden" / "admixtools2" / "v2.0.10" / "dg_6pop" / "f2_dir")
LEFT = ["French.DG", "Han.DG", "Papuan.DG"]
RIGHT = ["Mbuti.DG", "Yoruba.DG", "Chimp.REF"]


@pytest.fixture(scope="module")
def dg6():
    return pm.read_admixtools2_f2(DG6)


@pytest.fixture(scope="module")
def family6(tmp_path_factory) -> str:
    from test_skill import example_data

    d = tmp_path_factory.mktemp("family6")
    example_data(d)
    return str(d / "family6")


def popmix(*args: str) -> str:
    out = subprocess.run(["popmix", *args], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out.stdout


def test_alpha_sets_the_level_of_the_rank_tests(dg6):
    # rank 1 has p = 0.067: not rejected at 0.05, rejected at 0.1
    r = pm.qpwave(dg6, LEFT, RIGHT)
    assert r.rank_tests[1].p == pytest.approx(0.0670243, rel=1e-5)
    assert r.alpha == 0.05
    assert r.estimated_rank == 1
    r = pm.qpwave(dg6, LEFT, RIGHT, alpha=0.1)
    assert r.alpha == 0.1
    assert r.estimated_rank == 2


def test_the_result_makes_no_clade_claim(dg6):
    assert not hasattr(pm.qpwave(dg6, LEFT, RIGHT), "is_clade")


def test_cli_takes_alpha_and_reports_it(dg6, tmp_path):
    f2_dir = str(tmp_path / "f2")
    dg6.save(f2_dir)
    args = ["qpwave", "--f2-dir", f2_dir, "--left", ",".join(LEFT), "--right", ",".join(RIGHT), "--format", "json"]
    out = json.loads(popmix(*args, "--alpha", "0.1"))
    assert out["alpha"] == 0.1
    assert out["estimated_rank"] == 2
    assert "is_clade" not in out


def test_cli_calls_the_rank_a_lower_bound(dg6, tmp_path):
    f2_dir = str(tmp_path / "f2")
    dg6.save(f2_dir)
    text = popmix("qpwave", "--f2-dir", f2_dir, "--left", ",".join(LEFT), "--right", ",".join(RIGHT))
    assert "Rank tests (* rejected at alpha 0.05):" in text
    assert "Estimated rank: 1 (the lowest rank not rejected)" in text
    assert ("=> At least 2 streams of ancestry, as seen from these right populations: a lower bound, as other "
            "right populations or more data can show more streams (the rank can be at most 2 here).") in " ".join(
        text.split())


def test_cli_does_not_call_rank_0_a_clade(family6):
    # French_1 and French_2: rank 0 has p = 0.027, not rejected at 0.01
    text = popmix("qpwave", "-i", family6, "--left", "French_1,French_2", "--right", "Yoruba_1,Yoruba_2,Yoruba_3",
                  "--alpha", "0.01")
    assert "Estimated rank: 0 (the lowest rank not rejected)" in text
    assert ("=> Rank 0 not rejected: these right populations do not distinguish the left ones. That does not "
            "show that they form a clade: other right populations or more data can show more streams (the rank "
            "can be at most 1 here).") in " ".join(text.split())
