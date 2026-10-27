"""z is NaN, not 0 or infinite, when the standard error is 0 (G8, ADMIXTOOLS2-0008).

A statistic that is 0 by construction, such as f4(A, A; B, C), has SE 0 and
no z; the same goes for any estimate whose SE is 0.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
DG6 = str(ROOT / "tests" / "golden" / "admixtools2" / "v2.0.10" / "dg_6pop" / "f2_dir")


@pytest.fixture(scope="module")
def dg6():
    return pm.read_admixtools2_f2(DG6)


@pytest.mark.parametrize("pops", [
    ("Mbuti.DG", "Mbuti.DG", "Han.DG", "French.DG"),
    ("Mbuti.DG", "Han.DG", "French.DG", "French.DG"),
], ids=["pop1=pop2", "pop3=pop4"])
def test_f4_of_a_population_with_itself_has_no_z(dg6, pops):
    r = dg6.f4(*pops)
    assert r["est"] == 0.0
    assert r["se"] == 0.0
    assert math.isnan(r["z"])
    assert math.isnan(r["p"])


def test_qpdstat_cli_prints_no_z_without_an_se(dg6, tmp_path):
    f2_dir = str(tmp_path / "f2")
    dg6.save(f2_dir)
    out = subprocess.run(["popmix", "qpdstat", "--f2-dir", f2_dir, "--pop1", "Mbuti.DG", "--pop2", "Mbuti.DG",
                          "--pop3", "Han.DG", "--pop4", "French.DG", "--format", "json"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    row = json.loads(out.stdout)[0]
    assert row["est"] == 0.0
    assert row["se"] == 0.0
    assert row["z"] is None


def test_compare_fits_has_no_z_when_every_replicate_agrees():
    r = pm.compare_fits([2.0, 3.0, 4.0], [1.0, 2.0, 3.0])
    assert r["diff"] == 1.0
    assert r["se"] == 0.0
    assert math.isnan(r["z"])
    assert math.isnan(r["p"])
