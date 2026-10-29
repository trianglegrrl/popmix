"""The README is the short version; docs/guide.md has the detail (TASKS.md L7)."""

from pathlib import Path

ROOT = Path(__file__).parents[2]
README = (ROOT / "README.md").read_text(encoding="utf-8")
GUIDE_PATH = ROOT / "docs" / "guide.md"


def test_the_readme_is_brief():
    assert len(README.splitlines()) <= 250


def test_the_readme_links_the_guide():
    assert "(docs/guide.md)" in README


def test_the_guide_walks_through_the_commands():
    guide = GUIDE_PATH.read_text(encoding="utf-8")
    for command in ("info", "f2", "f3", "f4", "qpadm", "qpwave", "dates", "convertf"):
        assert f"### `popmix {command}`" in guide, command


def test_the_readme_uses_canadian_spelling():
    for american in ("licenses", "behavior", "color", "center"):
        assert american not in README.lower().replace("thirdparty", ""), american
