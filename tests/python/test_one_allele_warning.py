"""f2 and f3 say when a population's sampling bias cannot be corrected (G7).

A population of one pseudo-haploid sample has one allele per SNP, so the
bias p(1 - p)/(n - 1) of f2 with it, and of f3 with it as the target, cannot
be estimated: on f2 blocks they come out too large (admixtools2 computes
them the same way), where qp3pop on genotypes gives NaN with n = 0. popmix
warns on the f2 path; f4 and f3 with it as a source are not affected.
"""

from __future__ import annotations

import math
import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = ROOT / "tests" / "golden" / "dates" / "family_chr20_22"
POPS = ["Simulation", "French", "Yoruba"]


def message(statistic: str) -> str:
    return (f"Simulation has at most one allele at every SNP (as a single pseudo-haploid sample): its sampling bias "
            f"cannot be corrected, so {statistic} comes out too large (f4, and f3 with it as a source, are not "
            "affected)")


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


@pytest.fixture(scope="module")
def lone(tmp_path_factory) -> str:
    """The small DATES data with Simulation, one sample, made pseudo-haploid."""
    prefix = str(tmp_path_factory.mktemp("lone") / "lone")
    subprocess.run(["popmix", "convertf", "-i", str(SMALL), "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    lines = Path(f"{prefix}.geno").read_text(encoding="utf-8").splitlines()
    assert Path(f"{prefix}.ind").read_text(encoding="utf-8").split()[2] == "Simulation"
    Path(f"{prefix}.geno").write_text("".join(("2" if x[0] == "1" else x[0]) + x[1:] + "\n" for x in lines), encoding="utf-8")
    return prefix


@pytest.fixture(scope="module")
def f2(lone):
    return pm.extract_f2(lone, POPS)


def test_qp3pop_on_genotypes_has_no_value_for_such_a_target(lone):
    row = pm.qp3pop(lone, "Simulation", "French", "Yoruba").iloc[0]
    assert math.isnan(row["est"])
    assert row["n"] == 0


def test_f3_with_it_as_the_target_warns(f2):
    assert caught(lambda: f2.f3("Simulation", "French", "Yoruba")) == [message("f3 with it as the target")]


def test_f3_with_it_as_a_source_does_not_warn(f2):
    assert caught(lambda: f2.f3("French", "Simulation", "Yoruba")) == []


def test_f2_with_it_warns(f2):
    assert caught(lambda: f2.f2("Simulation", "French")) == [message("f2 with it")]


def test_f2_blocks_from_admixtools2_cannot_tell():
    dg6 = pm.read_admixtools2_f2(str(ROOT / "tests/golden/admixtools2/v2.0.10/dg_6pop/f2_dir"))
    assert caught(lambda: dg6.f3("Chimp.REF", "Mbuti.DG", "Han.DG")) == []


def test_the_cli_f3_warns(lone):
    out = subprocess.run(["popmix", "f3", "-i", lone, "--target", "Simulation", "--source1", "French",
                          "--source2", "Yoruba"], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert f"Warning: {message('f3 with it as the target')}" in out.stderr
