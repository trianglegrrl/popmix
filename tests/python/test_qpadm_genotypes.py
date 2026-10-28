"""qpAdm and qpWave with each f4 from the genotypes (TASKS.md H11).

``pm.qpadm(prefix, ...)`` is admixtools2's ``qpadm(prefix)``; the CLI takes
it as ``--f4-snps model`` (or ``each``, ``allsnps = TRUE``). The numbers are
checked against admixtools2 in tests/golden/test_qpadm_geno_admixtools2.py.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")
MODEL = ["--target", "Simulation", "--left", "French", "--right", "Yoruba,Simulation.1"]


def popmix(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["popmix", *args], capture_output=True, text=True)


@pytest.mark.parametrize(("flag", "allsnps"), [("model", False), ("each", True)])
def test_the_cli_fits_on_genotypes(flag, allsnps):
    r = popmix("qpadm", "-i", SMALL, *MODEL, "--f4-snps", flag, "--format", "json")
    assert r.returncode == 0, r.stderr
    want = pm.qpadm(SMALL, ["French"], ["Yoruba", "Simulation.1"], "Simulation", allsnps=allsnps)
    d = json.loads(r.stdout)
    assert d["chisq"] == pytest.approx(want.chisq, rel=1e-12)
    assert (d["n_blocks"], d["n_snps"]) == (want.n_blocks, want.n_snps)


def test_model_snps_include_monomorphic_ones():
    # f2 blocks keep only SNPs polymorphic across the populations
    f2 = pm.qpadm(pm.extract_f2(SMALL, ["Simulation", "French", "Yoruba", "Simulation.1"]), ["French"],
                  ["Yoruba", "Simulation.1"], "Simulation")
    geno = pm.qpadm(SMALL, ["French"], ["Yoruba", "Simulation.1"], "Simulation")
    assert geno.n_snps > f2.n_snps


def test_qpwave_on_genotypes():
    r = popmix("qpwave", "-i", SMALL, "--left", "Simulation,French", "--right", "Yoruba,Simulation.1",
               "--f4-snps", "model", "--format", "json")
    assert r.returncode == 0, r.stderr
    want = pm.qpwave(SMALL, ["Simulation", "French"], ["Yoruba", "Simulation.1"])
    assert json.loads(r.stdout)["rank_tests"][0]["chisq"] == pytest.approx(want.rank_tests[0].chisq, rel=1e-12)


def test_genotype_f4_needs_genotypes(tmp_path):
    pm.extract_f2(SMALL, ["Simulation", "French", "Yoruba", "Simulation.1"]).save(str(tmp_path / "f2"))
    r = popmix("qpadm", "--f2-dir", str(tmp_path / "f2"), *MODEL, "--f4-snps", "model")
    assert r.returncode == 2
    assert "give --input" in r.stderr


def test_f2_options_do_not_apply_to_genotype_f4():
    r = popmix("qpadm", "-i", SMALL, *MODEL, "--f4-snps", "model", "--maxmiss", "0.5")
    assert r.returncode == 2
    assert "--maxmiss apply to f2 blocks" in r.stderr
