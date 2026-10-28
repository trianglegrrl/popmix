"""Data paths in a reference, read on this machine.

A reference records the data it was made from by absolute path, under the
reference cache of the machine that made it
(`/home/a/.cache/popmix-references/aadr_v66p1/aadr_v66p1`). The goldens read
that data from this machine's cache, $POPMIX_REFERENCE_CACHE or
`~/.cache/popmix-references`, where `tests/golden/aadr/archive_aadr_v66p1.py
fetch` or the build scripts put it.
"""

from __future__ import annotations

import os
from pathlib import Path

MARKER = "/.cache/popmix-references/"


def cache() -> Path:
    return Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references"))


def local(path: str) -> str:
    """`path` with a reference-cache prefix moved to this machine's cache."""
    if MARKER not in path:
        return path
    return str(cache() / path.split(MARKER, 1)[1])
