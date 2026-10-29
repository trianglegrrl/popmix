"""Data paths in a reference, read on this machine.

A reference records the data it was made from by absolute path, under the
reference cache of the machine that made it
(`/home/a/.cache/popmix-references/aadr_v66p1/aadr_v66p1`). The goldens read
that data from this machine's cache, $POPMIX_REFERENCE_CACHE or
`~/.cache/popmix-references`, where `tests/golden/aadr/archive_aadr_v66p1.py
fetch` or the build scripts put it.

A golden skips when its data are missing. Under POPMIX_REQUIRE_GOLDEN=1 (CI),
data the repository ships must be there: `missing` raises instead, so a
golden that ships its data cannot skip unnoticed (TASKS.md K1).
"""

from __future__ import annotations

import os
from pathlib import Path

MARKER = "/.cache/popmix-references/"
REPO = Path(__file__).resolve().parents[2]


class ShippedDataMissing(RuntimeError):
    """Data the repository ships are missing, under POPMIX_REQUIRE_GOLDEN=1."""


def cache() -> Path:
    return Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references"))


def local(path: str) -> str:
    """`path` with a reference-cache prefix moved to this machine's cache."""
    if MARKER not in path:
        return path
    return str(cache() / path.split(MARKER, 1)[1])


def missing(path: str | os.PathLike) -> bool:
    """Whether `path` (relative to the repository root) does not exist.

    Raises
    ------
    ShippedDataMissing
        If it is missing, lies inside the repository and
        POPMIX_REQUIRE_GOLDEN=1 is set.
    """
    full = (REPO / path).resolve()
    if full.exists():
        return False
    if os.environ.get("POPMIX_REQUIRE_GOLDEN") == "1" and full.is_relative_to(REPO):
        raise ShippedDataMissing(f"{full}: shipped golden data missing (POPMIX_REQUIRE_GOLDEN=1)")
    return True
