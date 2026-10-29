"""Projected samples with missing data are said to be pulled in (G6).

Projection by dot product counts a missing genotype as 0, which pulls a
sample with missing data toward the origin, so popmix warns and points to
least squares.
"""

from __future__ import annotations

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
def gappy(tmp_path_factory) -> str:
    """The small data with 30% of NA00000's (Simulation's) genotypes missing."""
    prefix = str(tmp_path_factory.mktemp("gappy") / "g")
    subprocess.run(["popmix", "convertf", "-i", str(SMALL), "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    lines = Path(f"{prefix}.geno").read_text(encoding="utf-8").splitlines()
    Path(f"{prefix}.geno").write_text("".join(("9" if k % 10 < 3 else x[0]) + x[1:] + "\n"
                                              for k, x in enumerate(lines)), encoding="utf-8")
    return prefix


def test_dot_product_projection_with_missing_data_warns(gappy):
    messages = caught(lambda: pm.smartpca(gappy, ["French", "Yoruba"], ["Simulation"], n_components=2,
                                          outlier_iterations=0))
    assert len(messages) == 1
    assert messages[0].startswith("1 projected sample(s) miss more than 10% of the PCA's SNPs (up to 30%): ")


def test_least_squares_projection_does_not_warn(gappy):
    assert caught(lambda: pm.smartpca(gappy, ["French", "Yoruba"], ["Simulation"], n_components=2,
                                      outlier_iterations=0, lsqproject=True)) == []


def test_the_result_reports_missing_fractions(gappy):
    r = pm.smartpca(gappy, ["French", "Yoruba"], ["Simulation"], n_components=2, outlier_iterations=0,
                    lsqproject=True)
    sim = r.coords.index[r.coords.id == "NA00000"][0]
    assert r.missing_fraction[sim] == pytest.approx(0.3, abs=0.01)
