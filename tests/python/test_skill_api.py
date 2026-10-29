"""The skill's Python API reference keeps up with the package (TASKS.md L8).

``python-api.md`` names every member of ``popmix.__all__``, and every field,
property and method of each result dataclass in that result's own section,
so a name or field added, renamed or removed without its documentation fails
here. The repository's ``skills/popmix`` is the same skill, not a copy.
"""

from __future__ import annotations

import dataclasses
import re
import sys
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).resolve().parents[2]
SKILL = ROOT / "python" / "popmix" / "skill"
API = (SKILL / "python-api.md").read_text(encoding="utf-8")

RESULTS = sorted(n for n in pm.__all__ if dataclasses.is_dataclass(getattr(pm, n)))
# results documented inside another's section, or a section of the file
SECTION_OF = {
    "RankTest": "### `QpWaveResult`",
    "PopDropRow": "### `QpAdmResult`",
    "QpAdmAssessment": "### `QpAdmResult`",
    "LdBin": "### `DatesResult`",
    "DecayFit": "### `DatesResult`",
    "PulseCheck": "### `DatesResult`",
    "QpGraphEdge": "### `QpGraphResult`",
    "QpGraphF3": "### `QpGraphResult`",
    "QpGraphStart": "### `QpGraphResult`",
    "ScoredGraph": "## Admixture graphs",
}


def section(heading: str) -> str:
    """The text under ``heading`` up to the next heading of its level or above."""
    level = heading.split(" ", 1)[0]
    start = API.index(heading + "\n")
    stop = re.compile(rf"^#{{1,{len(level)}}} ", re.MULTILINE).search(API, start + len(heading) + 1)
    return API[start:stop.start() if stop else len(API)]


def members(cls: type) -> list[str]:
    """A result's fields, and its public properties and methods."""
    names = [f.name for f in dataclasses.fields(cls)]
    names += [k for k, v in vars(cls).items()
              if not k.startswith("_") and (isinstance(v, property) or callable(v)) and k not in names]
    return names


@pytest.mark.parametrize("name", sorted(pm.__all__))
def test_every_public_name_is_documented(name):
    assert f"`{name}`" in API


def test_no_table_row_for_a_name_that_does_not_exist():
    for row in re.findall(r"^\| (`[^|]+`) \|", API, re.MULTILINE):
        for name in re.findall(r"`([^`]+)`", row):
            assert name in pm.__all__, name


@pytest.mark.parametrize("name", RESULTS)
def test_every_result_field_is_documented_in_its_section(name):
    text = section(SECTION_OF.get(name, f"### `{name}`"))
    for member in members(getattr(pm, name)):
        assert re.search(rf"`{re.escape(member)}(\(|`)", text), f"{name}.{member}"


def test_skill_md_points_to_the_api_reference():
    assert "(python-api.md)" in (SKILL / "SKILL.md").read_text(encoding="utf-8")


@pytest.mark.skipif(sys.platform == "win32" and not (ROOT / "skills" / "popmix").is_symlink(),
                    reason="git checked the symlink out as a file (core.symlinks false)")
def test_the_repository_skill_is_the_package_skill():
    linked = ROOT / "skills" / "popmix"
    assert linked.is_symlink()
    assert linked.resolve() == SKILL.resolve()
