"""Too few jackknife blocks or chromosomes is said, not hidden (TASKS.md G4).

f2 blocks warn below 20 blocks; qpAdm and qpWave warn when their f4
statistics are too many for the blocks (chi-squared treats the jackknife
covariance as known); DATES and ROLLOFF warn when their jackknife runs over
fewer than 10 chromosomes.
"""

from __future__ import annotations

import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


def test_f2_blocks_warn_when_they_are_few():
    messages = caught(lambda: pm.extract_f2(SMALL, POPS, block_size=0.5))
    n = pm.extract_f2(SMALL, POPS, block_size=0.5).n_blocks
    assert n < 20
    assert messages == [f"only {n} jackknife blocks: standard errors, z and p-values from fewer than 20 are rough "
                        "(use more data or a smaller block size)"]


def test_enough_f2_blocks_do_not_warn():
    assert pm.extract_f2(SMALL, POPS).n_blocks >= 20
    assert caught(lambda: pm.extract_f2(SMALL, POPS)) == []


@pytest.fixture(scope="module")
def family6(tmp_path_factory) -> str:
    from test_skill import example_data

    d = tmp_path_factory.mktemp("family6")
    example_data(d)
    return str(d / "family6")


FAMILY6 = ["French_1", "French_2", "French_3", "Simulation.1", "Yoruba_1", "Yoruba_2", "Yoruba_3"]
RIGHT = ["French_2", "Yoruba_2", "French_3", "Yoruba_3"]


def test_qpadm_warns_when_its_f4_statistics_are_many_for_the_blocks(family6):
    # 2 sources and 4 right populations: 2 x 3 = 6 f4 statistics
    f2 = pm.extract_f2(family6, FAMILY6, block_size=0.1)
    g = f2.n_blocks
    messages = caught(lambda: pm.qpadm(f2, ["French_1", "Yoruba_1"], RIGHT, "Simulation.1"))
    assert messages == [f"{g} jackknife blocks for 6 f4 statistics: chi-squared treats their jackknife covariance "
                        f"as known, which inflates it about {(g - 1) / (g - 8):.2f}-fold here, so p is too small "
                        "(use more blocks or fewer populations)"]


def test_qpwave_with_few_statistics_for_its_blocks_does_not_warn(family6):
    # 2 left, 3 right: 2 f4 statistics, inflation (g - 1)/(g - 4) below 1.1
    f2 = pm.extract_f2(family6, FAMILY6)
    assert (f2.n_blocks - 1) / (f2.n_blocks - 4) < 1.1
    assert caught(lambda: pm.qpwave(f2, ["French_1", "French_2"], ["Yoruba_1", "Yoruba_2", "Yoruba_3"])) == []


def test_qpwave_warns_when_the_blocks_cannot_estimate_the_covariance(family6):
    # 4 left and 3 right populations: 3 x 2 = 6 f4 statistics from 7 blocks
    f2 = pm.extract_f2(family6, FAMILY6, block_size=0.5)
    g = f2.n_blocks
    assert g <= 8
    messages = caught(lambda: pm.qpwave(f2, ["Simulation.1", "French_1", "Yoruba_1", "French_2"],
                                        ["Yoruba_2", "French_3", "Yoruba_3"]))
    assert messages == [f"{g} jackknife blocks for 6 f4 statistics: too few to estimate their covariance, so "
                        "chi-squared and p mean nothing (use more blocks or fewer populations)"]


def test_dates_warns_about_three_chromosomes():
    messages = caught(lambda: pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"]))
    assert messages == ["DATES: the jackknife over 3 chromosomes gives a rough standard error (fewer than 10)"]


def test_dates_without_the_jackknife_does_not_warn():
    assert caught(lambda: pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], jackknife=False)) == []


def test_rolloff_warns_about_three_chromosomes():
    messages = caught(lambda: pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"]))
    assert messages == ["ROLLOFF: the jackknife over 3 chromosomes gives a rough standard error (fewer than 10)"]


def test_the_cli_prints_the_warning_on_stderr():
    out = subprocess.run(["popmix", "dates", "-i", SMALL, "--admixed", "Simulation.1", "--sources", "French,Yoruba"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "Warning: DATES: the jackknife over 3 chromosomes gives a rough standard error (fewer than 10)" in \
        out.stderr
