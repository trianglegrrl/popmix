"""Projected samples' shrinkage toward the origin, said and corrected (G6).

A PCA sample's coordinates include its own noise, which fitted the axes; a
sample projected from the same population has none of it and sits closer to
the origin. EIGENSOFT's documentation states the test: two halves of one
population, one in the PCA and one projected, should not look different.
``shrink=True`` divides projected coordinates by each significant
component's shrinkage factor. Projection by dot product also pulls samples
with missing genotypes toward the origin, so popmix warns and points to
least squares.
"""

from __future__ import annotations

import math
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
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
def halves(tmp_path_factory) -> str:
    """French and Yoruba, each split: half in the PCA, half (_B) projected."""
    prefix = tmp_path_factory.mktemp("halves") / "d"
    for ext in ("geno", "snp"):
        Path(f"{prefix}.{ext}").symlink_to(f"{SMALL}.{ext}")
    rows = [line.split() for line in Path(f"{SMALL}.ind").read_text().splitlines()]
    seen: dict[str, int] = {}
    out = []
    for i, sex, pop in rows:
        k = seen[pop] = seen.get(pop, 0) + 1
        out.append(f"{i} {sex} {pop}_B\n" if pop in ("French", "Yoruba") and k % 2 == 0 else f"{i} {sex} {pop}\n")
    Path(f"{prefix}.ind").write_text("".join(out))
    return str(prefix)


def run(prefix: str, **kw) -> pm.PcaResult:
    return pm.smartpca(prefix, ["French", "Yoruba"], ["French_B", "Yoruba_B"], n_components=4,
                       outlier_iterations=0, **kw)


def pc1_gap_ratio(r: pm.PcaResult) -> float:
    m = r.coords.groupby("pop")["PC1"].mean()
    return (m["French_B"] - m["Yoruba_B"]) / (m["French"] - m["Yoruba"])


def test_projected_halves_are_shrunk(halves):
    assert pc1_gap_ratio(run(halves)) < 0.9


def test_shrink_puts_the_projected_half_with_its_pca_half(halves):
    assert abs(pc1_gap_ratio(run(halves, shrink=True)) - 1) < 0.05


def test_shrink_corrects_only_significant_components(halves):
    r = run(halves, shrink=True)
    assert 0.5 < r.shrinkage[0] < 1
    assert all(math.isnan(x) for x in r.shrinkage[1:])
    plain = run(halves)
    pca = r.coords.role == "pca"
    pd.testing.assert_frame_equal(r.coords[pca], plain.coords[pca])
    proj = ~pca
    np.testing.assert_allclose(r.coords[proj]["PC1"], plain.coords[proj]["PC1"] / r.shrinkage[0], rtol=1e-12)
    np.testing.assert_array_equal(r.coords[proj]["PC2"], plain.coords[proj]["PC2"])


def test_without_shrink_there_are_no_factors(halves):
    assert len(run(halves).shrinkage) == 0


@pytest.fixture(scope="module")
def gappy(tmp_path_factory) -> str:
    """The small data with 30% of NA00000's (Simulation's) genotypes missing."""
    prefix = str(tmp_path_factory.mktemp("gappy") / "g")
    subprocess.run(["popmix", "convertf", "-i", str(SMALL), "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    lines = Path(f"{prefix}.geno").read_text().splitlines()
    Path(f"{prefix}.geno").write_text("".join(("9" if k % 10 < 3 else x[0]) + x[1:] + "\n"
                                              for k, x in enumerate(lines)))
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


def test_the_cli_takes_shrink(halves, tmp_path):
    out = subprocess.run(["popmix", "pca", "-i", halves, "--pca-pops", "French,Yoruba",
                          "--project-pops", "French_B,Yoruba_B", "--components", "4", "--outlier-iterations", "0",
                          "--shrink", "-o", str(tmp_path / "p")], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    coords = pd.read_csv(tmp_path / "p.coords.tsv", sep="\t")
    r = run(halves, shrink=True)
    np.testing.assert_allclose(coords["PC1"], r.coords["PC1"], rtol=1e-10)
