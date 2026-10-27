"""Population names that are unknown, repeated or in two roles are errors.

A typo used to pass silently: extraction gave a population with no data,
every statistic with it returned None, and qpAdm with the target among the
right populations reported a good fit (admixtools2 stops with "Duplicated
pops"). Each case now raises ValueError naming the populations.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

import popmix as pm

HERE = Path(__file__).parent
SMALL = str(HERE.parent / "golden" / "dates" / "family_chr20_22")
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]
DG6 = str(HERE.parent / "golden" / "admixtools2" / "v2.0.10" / "dg_6pop" / "f2_dir")


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(SMALL, POPS)


@pytest.fixture(scope="module")
def dg6():
    return pm.read_admixtools2_f2(DG6)


@pytest.fixture(scope="module")
def text(tmp_path_factory) -> str:
    out = tmp_path_factory.mktemp("text") / "small"
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", str(out), "-of", "eigenstrat"], check=True,
                   capture_output=True)
    return str(out)


@pytest.mark.parametrize("which", ["packed", "text"])
def test_extracting_an_unknown_population_is_an_error(text, which):
    prefix = SMALL if which == "packed" else text
    with pytest.raises(ValueError, match="NoSuchPop"):
        pm.extract_f2(prefix, ["French", "Yoruba", "NoSuchPop"])


@pytest.mark.parametrize("which", ["packed", "text"])
def test_extracting_a_population_twice_is_an_error(text, which):
    prefix = SMALL if which == "packed" else text
    with pytest.raises(ValueError, match="French"):
        pm.extract_f2(prefix, ["French", "French", "Yoruba"])


@pytest.mark.parametrize("stat, args", [
    ("f2", ("French", "Nope")),
    ("f3", ("Nope", "French", "Yoruba")),
    ("f4", ("French", "Yoruba", "Simulation", "Nope")),
    ("fst", ("Nope", "French")),
])
def test_statistics_of_an_unknown_population_are_an_error(f2, stat, args):
    with pytest.raises(ValueError, match="Nope"):
        getattr(f2, stat)(*args)


@pytest.mark.parametrize("target, sources, right, named", [
    ("French.DG", ["Han.DG", "Yoruba.DG"], ["Mbuti.DG", "French.DG", "Chimp.REF"], "French.DG"),
    ("French.DG", ["Han.DG", "Yoruba.DG"], ["Mbuti.DG", "Han.DG", "Chimp.REF"], "Han.DG"),
    ("French.DG", ["Han.DG", "Han.DG"], ["Mbuti.DG", "Papuan.DG", "Chimp.REF"], "Han.DG"),
    ("French.DG", ["French.DG", "Han.DG"], ["Mbuti.DG", "Papuan.DG", "Chimp.REF"], "French.DG"),
    ("French.DG", ["Han.DG", "Yoruba.DG"], ["Mbuti.DG", "Mbuti.DG", "Chimp.REF"], "Mbuti.DG"),
])
def test_qpadm_populations_must_be_distinct(dg6, target, sources, right, named):
    with pytest.raises(ValueError, match=named):
        pm.qpadm(dg6, sources, right, target)


def test_qpadm_with_an_unknown_population_is_an_error(dg6):
    with pytest.raises(ValueError, match="Nope"):
        pm.qpadm(dg6, ["Han.DG", "Nope"], ["Mbuti.DG", "Papuan.DG", "Chimp.REF"], "French.DG")


def test_qpwave_left_and_right_must_be_distinct(dg6):
    with pytest.raises(ValueError, match="Han.DG"):
        pm.qpwave(dg6, ["French.DG", "Han.DG"], ["Han.DG", "Mbuti.DG", "Chimp.REF"])


def test_qpadm_rotation_with_an_unknown_target_is_an_error(dg6):
    with pytest.raises(ValueError, match="TypoTarget"):
        pm.qpadm_rotate(dg6, leftright=["Han.DG", "Yoruba.DG", "Mbuti.DG"], target="TypoTarget",
                        rightfix=["Papuan.DG", "Chimp.REF"])


def test_cli_qpadm_rejects_the_target_among_the_right_populations(dg6, tmp_path):
    d = tmp_path / "f2"
    dg6.save(str(d))
    out = subprocess.run(["popmix", "qpadm", "--f2-dir", str(d), "--target", "French.DG",
                          "--left", "Han.DG,Yoruba.DG", "--right", "Mbuti.DG,French.DG,Chimp.REF"],
                         capture_output=True, text=True)
    assert out.returncode != 0
    assert "French.DG" in out.stderr


def test_allele_frequencies_of_an_unknown_population_are_an_error():
    geno = pm.read_genotypes(SMALL)
    with pytest.raises(ValueError, match="Nope"):
        geno.allele_frequencies(populations=["French", "Nope"])
    af = geno.allele_frequencies(populations=["French", "Yoruba"])
    with pytest.raises(ValueError, match="Nope"):
        af.filter_populations(["Nope", "French"])
