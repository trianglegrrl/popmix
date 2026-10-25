"""The version and project URLs agree everywhere a release reads them.

PyPI reads pyproject.toml, the compiled module reads the Cargo workspace
version, and bioconda reads recipes/meta.yaml.
"""

import re
from pathlib import Path

import popmix as pm

ROOT = Path(__file__).resolve().parents[2]
REPOSITORY = "https://github.com/trianglegrrl/popmix"


def _match(text: str, pattern: str) -> str:
    found = re.search(pattern, text, re.MULTILINE)
    assert found, f"{pattern!r} not found"
    return found.group(1)


def _read(path: str) -> str:
    return (ROOT / path).read_text()


def test_versions_agree() -> None:
    pyproject = _match(_read("pyproject.toml"), r'^version = "([^"]+)"')
    cargo = _match(_read("Cargo.toml"), r'^version = "([^"]+)"')
    recipe = _match(_read("recipes/meta.yaml"), r'set version = "([^"]+)"')
    assert pyproject == cargo == recipe == pm.__version__


def test_project_urls_point_at_the_repository() -> None:
    urls = _read("pyproject.toml").split("[project.urls]")[1].split("\n[")[0]
    for key in ("Homepage", "Repository"):
        assert _match(urls, rf'^{key} = "([^"]+)"') == REPOSITORY
    recipe = _read("recipes/meta.yaml")
    for key in ("home", "dev_url"):
        assert _match(recipe, rf"^  {key}: (\S+)") == REPOSITORY

