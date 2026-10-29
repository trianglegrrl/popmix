"""Holdout on AADR v66.p1: 1000 Genomes populations, ancient groups by least squares.

The oracle is the test EIGENSOFT's documentation gives for projection: split
each PCA population by hand, put one half in the PCA and project the other.
``holdout=True`` must place every sample exactly as that run does, and its
measured shrinkage must be the slope of the projected halves' population
means on the PCA halves'. On these populations projected samples are pulled
in by a few percent on PC1 and PC2 and by a quarter to three quarters on
PC3-PC6, which is why popmix measures it rather than leaving it unsaid.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

AADR = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1" \
    / "aadr_v66p1"
KGP = ["CEU.DG", "TSI.DG", "FIN.DG", "YRI.DG", "CHB.DG", "GIH.DG"]
ANCIENT = ["Germany_HalberstadtSonntagsfeld_EN_LBK.AG", "Czechia_EBA_CordedWare.AG", "England_MIA.AG"]
PCS = [f"PC{j}" for j in range(1, 7)]

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not Path(f"{AADR}.geno").exists(), reason=f"AADR v66.p1 build not available: {AADR}"),
]


def run(prefix: str, pca: list[str], project: list[str], **kw) -> pm.PcaResult:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = pm.smartpca(prefix, pca, project, n_components=6, outlier_iterations=0, lsqproject=True, **kw)
    r.warnings = [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]
    return r


@pytest.fixture(scope="module")
def held() -> pm.PcaResult:
    return run(str(AADR), KGP, ANCIENT, holdout=True)


@pytest.fixture(scope="module")
def by_hand(tmp_path_factory) -> pm.PcaResult:
    """Every second sample of each 1000 Genomes population relabelled _B."""
    prefix = tmp_path_factory.mktemp("halves") / "d"
    for ext in ("geno", "snp"):
        Path(f"{prefix}.{ext}").symlink_to(f"{AADR}.{ext}")
    seen: dict[str, int] = {}
    out = []
    for i, sex, pop in (line.split() for line in Path(f"{AADR}.ind").read_text(encoding="utf-8").splitlines()):
        if pop in KGP:
            k = seen[pop] = seen.get(pop, 0) + 1
            pop = f"{pop}_B" if k % 2 == 0 else pop
        out.append(f"{i} {sex} {pop}\n")
    Path(f"{prefix}.ind").write_text("".join(out), encoding="utf-8")
    return run(str(prefix), KGP, [f"{p}_B" for p in KGP] + ANCIENT)


def test_held_out_samples_are_placed_as_projecting_them(held, by_hand):
    assert held.coords.id.tolist() == by_hand.coords.id.tolist()
    np.testing.assert_allclose(held.coords[PCS], by_hand.coords[PCS], rtol=1e-10, atol=1e-13)


def test_shrinkage_is_the_slope_of_the_projected_halves(held, by_hand):
    m = by_hand.coords.groupby("pop")[PCS].mean()
    a = m.loc[KGP].to_numpy()
    b = m.loc[[f"{p}_B" for p in KGP]].to_numpy()
    a, b = a - a.mean(0), b - b.mean(0)
    np.testing.assert_allclose(held.shrinkage, (a * b).sum(0) / (a * a).sum(0), rtol=1e-10)


def test_projected_samples_are_pulled_in_on_the_later_components(held):
    assert held.shrinkage[0] > 0.95
    assert held.shrinkage[5] < 0.5


def test_the_warning_names_every_component_pulled_in_by_more_than_a_tenth(held):
    low = [f"{s:.0%} (PC{j})" for j, s in enumerate(held.shrinkage, 1) if s < 0.9]
    assert len(low) >= 3
    assert held.warnings == ["projected samples sit closer to the origin than the PCA samples: held-out reference "
                             "samples reach " + ", ".join(low) + " of their populations' spread; compare "
                             "projected samples with the held-out ones (role \"holdout\"), not with the PCA "
                             "samples"]
