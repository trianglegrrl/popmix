"""No individual or population AADR v66.p1 removed appears in our data (O6).

v66.p1 withdrew the 161 present-day genomes of Jacobs et al. 2019 and the 27
populations made only of them; their data must not be used. Every token of
every tracked text file, and of the sample files in the reference cache, is
hashed and checked against tests/data_ethics/removed_v66p1.json
(tests/data_ethics/generate_removed.py), which stores only hashes.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
REMOVED = set(json.loads((ROOT / "tests/data_ethics/removed_v66p1.json").read_text())["sha256"])
CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references"))
# sample and population names: letters, digits, _ . - ( ); ':' joins a
# population and a sample ID in some .ind files
TOKEN = re.compile(r"[A-Za-z0-9_.()\-]+")
SAMPLE_FILES = {".ind", ".fam", ".pedind", ".json", ".txt", ".tsv", ".csv", ".ped"}


def _hits(text: str) -> set[str]:
    found = set()
    for token in set(TOKEN.findall(text)):
        for candidate in (token, token.rstrip(".")):
            if hashlib.sha256(candidate.encode()).hexdigest() in REMOVED:
                found.add(candidate)
    return found


def _text(path: Path) -> str | None:
    try:
        if path.suffix == ".ped":
            # a .ped line is family, sample, ... then every genotype: only the
            # names are needed
            with path.open() as f:
                return "\n".join(" ".join(line.split(maxsplit=2)[:2]) for line in f)
        return path.read_text()
    except (UnicodeDecodeError, OSError):
        return None


def test_the_hash_list_is_complete():
    assert len(REMOVED) == 376


def test_a_removed_name_would_be_found():
    """The check works: a name whose hash is listed is caught."""
    fake = "Pop_Example.DG"
    REMOVED.add(hashlib.sha256(fake.encode()).hexdigest())
    try:
        assert _hits(f"S1 M {fake}\nS2:{fake} F x\n") == {fake}
    finally:
        REMOVED.discard(hashlib.sha256(fake.encode()).hexdigest())


@pytest.mark.skipif(shutil.which("git") is None, reason="git not available")
def test_not_in_the_repository():
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True,
                           check=True).stdout.split("\0")
    found = {}
    for name in filter(None, files):
        text = _text(ROOT / name)
        if text and (hits := _hits(text)):
            found[name] = hits
    assert not found, f"removed v66.p1 data found: {sorted(found)}"


def test_not_in_the_reference_cache():
    if not CACHE.exists():
        pytest.skip("no reference cache")
    found = {}
    for path in CACHE.rglob("*"):
        if path.suffix in SAMPLE_FILES and path.is_file() and path.stat().st_size < 1 << 30:
            text = _text(path)
            if text and (hits := _hits(text)):
                found[str(path)] = hits
    assert not found, f"removed v66.p1 data found: {sorted(found)}"


# reads v62.0's annotation only to list the names v66.p1 removed
READS_V62_TO_HASH = {"tests/data_ethics/generate_removed.py"}


@pytest.mark.skipif(shutil.which("git") is None, reason="git not available")
def test_no_example_points_at_a_release_with_the_withdrawn_genomes():
    """Docs and examples name v66.p1 data, never a v62.0 or v66.0 path."""
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True,
                           check=True).stdout.split("\0")
    path = re.compile(r"v6[26]\.0_1240[kK]")
    found = sorted(name for name in filter(None, files)
                   if name not in READS_V62_TO_HASH and (text := _text(ROOT / name)) and path.search(text))
    assert found == []
