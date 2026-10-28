"""References name their data under the reference cache of the machine that
made them; the goldens read it from this machine's cache (TASKS.md K2)."""

from pathlib import Path

from .reference_paths import local


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
                 for n, line in enumerate(p.read_text().splitlines(), 1) if bare.search(line)]
    assert offenders == []
