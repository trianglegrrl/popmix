"""Kinship says when its estimate cannot be trusted (TASKS.md G11).

KING-robust needs heterozygous calls: a pseudo-haploid sample has none and
its kinship is minus infinity (as PLINK 2), so popmix warns and points to
READ. READ's fixed thresholds take its normalization (the median or mean P0)
to be that of unrelated pairs: with few pairs it warns.
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


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


@pytest.fixture(scope="module")
def lone(tmp_path_factory) -> str:
    """The small DATES data with NA00000 (Simulation) made pseudo-haploid."""
    prefix = str(tmp_path_factory.mktemp("lone") / "lone")
    subprocess.run(["popmix", "convertf", "-i", str(SMALL), "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    lines = Path(f"{prefix}.geno").read_text(encoding="utf-8").splitlines()
    Path(f"{prefix}.geno").write_text("".join(("2" if x[0] == "1" else x[0]) + x[1:] + "\n" for x in lines), encoding="utf-8")
    return prefix


def test_king_warns_about_samples_without_heterozygotes(lone):
    messages = caught(lambda: pm.king(lone, ["Simulation", "Simulation.1"]))
    assert messages == ["KING-robust needs heterozygous calls: 1 sample(s) have none (pseudo-haploid?), so their "
                        "kinship is minus infinity, as in PLINK 2 (NA00000); use READ (--method read) for "
                        "pseudo-haploid data"]
    table = pm.king(lone, ["Simulation", "Simulation.1"])
    assert all(math.isinf(k) for k in table[table.id2 == "NA00000"].kinship)


def test_king_on_diploid_samples_does_not_warn():
    assert caught(lambda: pm.king(str(SMALL), ["Simulation.1"])) == []


def test_read_with_few_pairs_warns():
    # 4 samples: 6 pairs
    messages = caught(lambda: pm.read_kinship(str(SMALL), ["Simulation.1"]))
    assert messages == ["READ: normalized by 6 pairs (fewer than 10); its thresholds take the norm to be unrelated "
                        "pairs, so with this few, or with many relatives, the degrees are unreliable (give a norm "
                        "value from a larger set)"]


def test_read_with_a_given_norm_value_does_not_warn():
    assert caught(lambda: pm.read_kinship(str(SMALL), ["Simulation.1"], norm_value=0.25)) == []


def test_the_cli_says_so(lone):
    out = subprocess.run(["popmix", "kinship", "-i", lone, "--populations", "Simulation,Simulation.1"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "Warning: KING-robust needs heterozygous calls: 1 sample(s) have none" in out.stderr
