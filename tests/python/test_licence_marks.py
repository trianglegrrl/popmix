"""Test tooling that calls admixtools2's internals is GPL-3, marked as such;
nothing that ships does (TASKS.md C11)."""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPDX = "# SPDX-License-Identifier: GPL-3.0-or-later"
INTERNAL = "admixtools:::"
CALLERS = sorted(p for p in (ROOT / "tests").rglob("*.R") if INTERNAL in p.read_text(encoding="utf-8"))


def test_some_generator_calls_admixtools2_internals() -> None:
    assert CALLERS != []


@pytest.mark.parametrize("path", CALLERS, ids=lambda p: p.name)
def test_a_script_calling_admixtools2_internals_is_marked_gpl(path: Path) -> None:
    assert SPDX in path.read_text(encoding="utf-8").splitlines()[:3]


@pytest.mark.parametrize("path", CALLERS, ids=lambda p: p.name)
def test_the_gpl_text_is_beside_it(path: Path) -> None:
    assert (path.parent / "COPYING").read_text(encoding="utf-8").startswith(
        "                    GNU GENERAL PUBLIC LICENSE\n                       Version 3, 29 June 2007")


def test_no_shipped_code_calls_admixtools2_internals() -> None:
    shipped = [p for d in ("python", "crates") for p in (ROOT / d).rglob("*")
               if p.suffix in (".py", ".rs", ".R", ".md") and INTERNAL in p.read_text(encoding="utf-8")]
    assert shipped == []
