"""DATES's date in years is honest about what it is (TASKS.md G3).

Generations are counted before the admixed samples lived: years before
present need their age. The generation time is an assumption with its own
uncertainty, and R² measures how well one exponential fits, not whether
there was a single pulse.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")
DATES = ["dates", "-i", SMALL, "--admixed", "Simulation.1", "--sources", "French,Yoruba"]


@pytest.fixture(scope="module")
def plain():
    return pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"])


def popmix(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["popmix", *args], capture_output=True, text=True)


def test_years_add_the_sample_age_and_every_uncertainty(plain):
    r = pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], generation_time=25.0, generation_time_se=3.0,
                 sample_age=4500.0, sample_age_se=50.0)
    g, g_se = plain.generations, plain.generations_se
    assert r.generations == pytest.approx(g, rel=1e-12)
    assert r.years_bp == pytest.approx(4500.0 + 25.0 * g, rel=1e-12)
    expected_se = math.sqrt(50.0 ** 2 + g ** 2 * 3.0 ** 2 + 25.0 ** 2 * g_se ** 2 + g_se ** 2 * 3.0 ** 2)
    assert r.years_bp_se == pytest.approx(expected_se, rel=1e-12)


def test_a_negative_sample_age_is_an_error():
    with pytest.raises(pm.InputError, match="sample_age"):
        pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], sample_age=-1.0, jackknife=False)


def test_cli_json_carries_the_assumptions(plain):
    out = popmix(*DATES, "--generation-time", "25", "--generation-time-se", "3", "--sample-age", "4500",
                 "--sample-age-se", "50", "--format", "json")
    assert out.returncode == 0, out.stderr
    d = json.loads(out.stdout)
    assert d["years_bp"] == pytest.approx(4500.0 + 25.0 * plain.generations, rel=1e-12)
    assert (d["generation_time"], d["generation_time_se"], d["sample_age"], d["sample_age_se"]) == \
        (25.0, 3.0, 4500.0, 50.0)


def test_cli_without_a_sample_age_dates_before_the_samples(plain):
    out = popmix(*DATES)
    assert out.returncode == 0, out.stderr
    years, se = plain.years_bp, plain.years_bp_se
    assert (f"Date: {years:.0f} ± {se:.0f} years before the admixed samples lived (29 years/gen, taken as "
            "exact; give --sample-age for years BP)") in " ".join(out.stdout.split())
    assert "Estimated years BP" not in out.stdout


def test_cli_with_a_sample_age_dates_before_present():
    r = pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], sample_age=4500.0, generation_time_se=3.0)
    out = popmix(*DATES, "--sample-age", "4500", "--generation-time-se", "3")
    assert out.returncode == 0, out.stderr
    assert (f"Date: {r.years_bp:.0f} ± {r.years_bp_se:.0f} years BP (sample age 4500 ± 0, "
            "29 ± 3 years/gen)") in " ".join(out.stdout.split())


def test_cli_does_not_present_r_squared_as_a_single_pulse():
    out = popmix(*DATES, "--no-jackknife")
    assert out.returncode == 0, out.stderr
    line = next(line for line in out.stdout.splitlines() if line.startswith("Fit R²"))
    assert line.endswith("(how well one exponential fits; it does not show a single pulse)")
