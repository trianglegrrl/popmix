"""References name their data under the reference cache of the machine that
made them; the goldens read it from this machine's cache (TASKS.md K2)."""

from pathlib import Path

import pytest

from .reference_paths import REPO, ShippedDataMissing, local, missing


def test_a_cache_path_moves_to_this_machines_cache(monkeypatch, tmp_path):
    monkeypatch.setenv("POPMIX_REFERENCE_CACHE", str(tmp_path))
    assert local("/home/a/.cache/popmix-references/aadr_v66p1/aadr_v66p1") == str(tmp_path / "aadr_v66p1/aadr_v66p1")


def test_the_default_cache_is_in_the_home_directory(monkeypatch):
    monkeypatch.delenv("POPMIX_REFERENCE_CACHE", raising=False)
    expected = Path.home() / ".cache/popmix-references/vcf_ped/family"
    assert local("/home/someone/.cache/popmix-references/vcf_ped/family") == str(expected)


def test_a_path_outside_the_cache_is_kept(monkeypatch, tmp_path):
    monkeypatch.setenv("POPMIX_REFERENCE_CACHE", str(tmp_path))
    assert local("tests/data/family/family") == "tests/data/family/family"


def test_every_golden_reads_a_reference_prefix_through_local():
    import re

    bare = re.compile(r"""(?<!local\()\b(?:f2_)?ref\[["']prefix["']\]""")
    offenders = [f"{p.name}:{n}" for p in sorted(Path(__file__).parent.glob("test_*.py"))
                 for n, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1) if bare.search(line)]
    assert offenders == []


def test_data_outside_the_repository_may_be_missing_even_when_goldens_are_required(monkeypatch, tmp_path):
    monkeypatch.setenv("POPMIX_REQUIRE_GOLDEN", "1")
    assert missing(tmp_path / "aadr_v66p1.geno") is True


def test_missing_shipped_data_lets_a_golden_skip_by_default(monkeypatch):
    monkeypatch.delenv("POPMIX_REQUIRE_GOLDEN", raising=False)
    assert missing(REPO / "tests/golden/no_such_reference.json") is True


def test_missing_shipped_data_fails_when_goldens_are_required(monkeypatch):
    monkeypatch.setenv("POPMIX_REQUIRE_GOLDEN", "1")
    with pytest.raises(ShippedDataMissing, match="no_such_reference.json"):
        missing(REPO / "tests/golden/no_such_reference.json")


def test_a_relative_path_is_read_from_the_repository(monkeypatch):
    monkeypatch.setenv("POPMIX_REQUIRE_GOLDEN", "1")
    with pytest.raises(ShippedDataMissing):
        missing("tests/golden/dates/no_such_prefix.geno")


def test_shipped_data_that_is_there_is_not_missing(monkeypatch):
    monkeypatch.setenv("POPMIX_REQUIRE_GOLDEN", "1")
    assert missing(REPO / "tests/golden/dates/family_chr20_22.geno") is False
