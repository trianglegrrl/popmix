"""The warnings channel (TASKS.md G13): a warning from the Rust core is a
popmix.PopmixWarning in Python and one ``Warning: ...`` line from the CLI."""

import subprocess
import sys
import warnings

import pytest

import popmix as pm

MULTIALLELIC = "tests/fixtures/small_vcf_multiallelic.vcf"


def test_popmix_warnings_are_user_warnings_of_their_own_class():
    assert issubclass(pm.PopmixWarning, UserWarning)
    assert "PopmixWarning" in pm.__all__
    assert pm.PopmixWarning.__module__ == "popmix.errors"


def test_a_warning_from_the_core_is_a_popmix_warning():
    with pytest.warns(pm.PopmixWarning, match=r"skipped 4 records that are not biallelic SNPs"):
        pm.read_vcf(MULTIALLELIC)


def test_popmix_warnings_can_be_silenced_or_made_errors_by_class():
    with warnings.catch_warnings():
        warnings.simplefilter("error", pm.PopmixWarning)
        with pytest.raises(pm.PopmixWarning):
            pm.read_vcf(MULTIALLELIC)
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("ignore", pm.PopmixWarning)
        pm.read_vcf(MULTIALLELIC)
    assert seen == []


def _cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "popmix.cli", *args], capture_output=True, text=True)


def test_the_cli_prints_a_warning_as_one_line():
    run = _cli("info", MULTIALLELIC)
    assert run.returncode == 0
    lines = run.stderr.splitlines()
    assert lines == [f"Warning: {MULTIALLELIC}: skipped 4 records that are not biallelic SNPs "
                     "(1 multiallelic, 3 indels or symbolic alleles)"]
    assert "Format: vcf" in run.stdout
