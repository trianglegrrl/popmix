"""Holding out reference samples, to compare projected samples like with like (G6).

A PCA sample's coordinates include its own noise, which fitted the axes; a
projected sample has none of it and sits closer to the origin. EIGENSOFT's
documentation gives the test: two halves of one population, one in the PCA
and one projected, should not look different, and they do. ``holdout=True``
projects every second sample of each PCA population, so projected samples
can be compared with reference samples placed the same way, and measures
each component's shrinkage from the held-out halves.
"""

from __future__ import annotations

import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = ROOT / "tests" / "golden" / "dates" / "family_chr20_22"
PCS = ["PC1", "PC2", "PC3", "PC4"]


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


def every_second(pops: tuple[str, ...]) -> list[str]:
    """Ids of the 2nd, 4th, ... sample of each population, in .ind order."""
    seen: dict[str, int] = {}
    out = []
    for i, _, pop in (line.split() for line in Path(f"{SMALL}.ind").read_text(encoding="utf-8").splitlines()):
        k = seen[pop] = seen.get(pop, 0) + 1
        if pop in pops and k % 2 == 0:
            out.append(i)
    return out


@pytest.fixture(scope="module")
def halves(tmp_path_factory) -> str:
    """French and Yoruba relabelled by hand: every second sample _B."""
    prefix = tmp_path_factory.mktemp("halves") / "d"
    for ext in ("geno", "snp"):
        Path(f"{prefix}.{ext}").symlink_to(f"{SMALL}.{ext}")
    held = set(every_second(("French", "Yoruba")))
    rows = [line.split() for line in Path(f"{SMALL}.ind").read_text(encoding="utf-8").splitlines()]
    Path(f"{prefix}.ind").write_text("".join(f"{i} {sex} {pop}_B\n" if i in held else f"{i} {sex} {pop}\n"
                                             for i, sex, pop in rows), encoding="utf-8")
    return str(prefix)


def run(**kw) -> pm.PcaResult:
    return pm.smartpca(SMALL, ["French", "Yoruba"], ["Simulation.1"], n_components=4, outlier_iterations=0, **kw)


@pytest.fixture(scope="module")
def held() -> pm.PcaResult:
    return run(holdout=True)


@pytest.fixture(scope="module")
def by_hand(halves) -> pm.PcaResult:
    return pm.smartpca(halves, ["French", "Yoruba"], ["French_B", "Yoruba_B", "Simulation.1"], n_components=4,
                       outlier_iterations=0)


def test_every_second_reference_sample_is_held_out(held):
    assert held.coords.id[held.coords.role == "holdout"].tolist() == every_second(("French", "Yoruba"))
    assert held.coords.id[held.coords.role == "project"].tolist() == \
        held.coords.id[held.coords["pop"] == "Simulation.1"].tolist()


def test_held_out_samples_are_placed_as_projecting_them(held, by_hand):
    np.testing.assert_allclose(held.coords[PCS], by_hand.coords[PCS], rtol=1e-10, atol=1e-13)
    assert held.n_pca_samples == by_hand.n_pca_samples


def test_shrinkage_is_the_held_out_gap_over_the_pca_gap(held, by_hand):
    m = by_hand.coords.groupby("pop")[PCS].mean()
    gap = (m.loc["French_B"] - m.loc["Yoruba_B"]) / (m.loc["French"] - m.loc["Yoruba"])
    np.testing.assert_allclose(held.shrinkage, gap.to_numpy(), rtol=1e-10)


def test_projected_halves_are_shrunk_on_pc1(held):
    assert held.shrinkage[0] < 0.9


def test_strong_shrinkage_warns():
    messages = caught(lambda: run(holdout=True))
    r = run(holdout=True)
    low = [(j, s) for j, s in enumerate(r.shrinkage, 1) if s < 0.9]
    assert low[0][0] == 1
    assert messages == ["projected samples sit closer to the origin than the PCA samples: held-out reference "
                        "samples reach " + ", ".join(f"{s:.0%} (PC{j})" for j, s in low)
                        + " of their populations' spread; compare projected samples with the held-out ones "
                        "(role \"holdout\"), not with the PCA samples"]


def test_without_holdout_nothing_is_measured():
    r = run()
    assert len(r.shrinkage) == 0
    assert "holdout" not in set(r.coords.role)


def test_holdout_needs_two_populations_to_measure():
    with pytest.raises(pm.InputError, match="holdout"):
        pm.smartpca(SMALL, ["French"], ["Yoruba"], n_components=2, outlier_iterations=0, holdout=True)


def test_the_cli_marks_and_reports_the_held_out_samples(held, tmp_path):
    out = subprocess.run(["popmix", "pca", "-i", str(SMALL), "--pca-pops", "French,Yoruba",
                          "--project-pops", "Simulation.1", "--components", "4", "--outlier-iterations", "0",
                          "--holdout", "-o", str(tmp_path / "p")], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    coords = pd.read_csv(tmp_path / "p.coords.tsv", sep="\t")
    assert coords.role.tolist() == held.coords.role.tolist()
    assert "Held-out reference samples reach " + ", ".join(f"{s:.0%} (PC{j})" for j, s in
                                                           enumerate(held.shrinkage, 1)) \
        + " of their PCA halves' spread" in out.stdout


def test_the_plot_draws_held_out_samples_apart(held):
    pytest.importorskip("matplotlib").use("Agg")
    from popmix import plot

    ax = plot.pca(held)
    # one collection per population and role: French and Yoruba twice, Simulation.1 once
    assert len(ax.collections) == 5
