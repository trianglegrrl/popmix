"""Every set of references records the environment that made it (TASKS.md K2):
`environment.txt` in each tool-version directory, written by
`record_environment.py`, names the tool's version and the R session or the
binary's checksum."""

from pathlib import Path

import pytest

from .record_environment import TOOLS

GOLDEN = Path(__file__).parent


def test_every_versioned_reference_directory_has_a_tool_entry():
    versioned = {str(p.relative_to(GOLDEN)) for p in GOLDEN.glob("*/*")
                 if p.is_dir() and any(p.rglob("*.json")) and p.name[:1] in "vr0123456789"
                 and p.name != "vcf_out"}
    assert versioned == set(TOOLS)


@pytest.mark.parametrize("directory", sorted(TOOLS))
def test_the_environment_is_recorded(directory):
    assert (GOLDEN / directory / "environment.txt").exists()


@pytest.mark.parametrize("directory", sorted(TOOLS))
def test_the_record_names_the_version_in_the_directory_name(directory):
    version = Path(directory).name.removeprefix("v").removeprefix("r-")
    assert version in (GOLDEN / directory / "environment.txt").read_text()


@pytest.mark.parametrize("directory", [d for d, t in TOOLS.items() if t.get("r")])
def test_an_r_reference_records_its_session(directory):
    assert "sessionInfo()" in (GOLDEN / directory / "environment.txt").read_text()
