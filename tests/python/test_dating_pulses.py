"""Is one admixture pulse enough? DATES and ROLLOFF check (TASKS.md H8).

Both fit two exponentials beside the single one, with jackknife SEs over
chromosomes, and support two pulses only if both amplitudes and the
dates' difference clear 2 SE. The small data's admixed sample was
simulated with one pulse, and is not taken for two; the fit's recovery of
two pulses and the decision on curves with known answers are unit-tested
in crates/popmix-core/src/dates/two_pulse.rs.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")


@pytest.fixture(scope="module")
def dated() -> pm.DatesResult:
    return pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"])


def test_dates_checks_for_a_second_pulse(dated):
    c = dated.pulse_check
    assert isinstance(c, pm.PulseCheck)
    assert c.generations[0] < c.generations[1]
    assert c.multiple_pulses is False


def test_the_check_needs_the_jackknife():
    assert pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], jackknife=False).pulse_check is None


def test_rolloff_checks_too():
    r = pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"])
    assert isinstance(r.pulse_check, pm.PulseCheck)
    assert r.pulse_check.multiple_pulses is False


def popmix(*args: str) -> subprocess.CompletedProcess:
    r = subprocess.run(["popmix", *args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


def test_the_cli_says_one_pulse_is_not_rejected(dated):
    out = " ".join(popmix("dates", "-i", SMALL, "--admixed", "Simulation.1", "--sources", "French,Yoruba").stdout.split())
    assert "One pulse is not rejected: a second exponential is not supported" in out


def test_the_json_carries_the_check(dated):
    d = json.loads(popmix("dates", "-i", SMALL, "--admixed", "Simulation.1", "--sources", "French,Yoruba",
                          "--format", "json").stdout)
    assert d["pulse_check"]["multiple_pulses"] is False
    assert d["pulse_check"]["generations"] == pytest.approx(list(dated.pulse_check.generations), rel=1e-12)
