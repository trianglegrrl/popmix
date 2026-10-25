"""The version and project URLs agree everywhere a release reads them.

PyPI reads pyproject.toml, the compiled module reads the Cargo workspace
version, and bioconda reads recipes/meta.yaml.
"""

import re
import shutil
import subprocess
import tarfile
from pathlib import Path

import pytest

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



@pytest.mark.slow
def test_sdist_carries_what_downstream_builds_need(tmp_path: Path) -> None:
    """bioconda builds from the sdist: it needs the license, the lockfile
    and every crate."""
    maturin = shutil.which("maturin")
    if maturin is None:
        pytest.skip("maturin not installed")
    subprocess.run(
        [maturin, "sdist", "-o", str(tmp_path)],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    (sdist,) = tmp_path.glob("popmix-*.tar.gz")
    with tarfile.open(sdist) as tar:
        names = {n.split("/", 1)[1] for n in tar.getnames() if "/" in n}
    for needed in (
        "LICENSE",
        "README.md",
        "Cargo.lock",
        "pyproject.toml",
        "crates/popmix-core/Cargo.toml",
        "crates/popmix-io/Cargo.toml",
        "crates/popmix-python/Cargo.toml",
        "python/popmix/__init__.py",
    ):
        assert needed in names, f"{needed} missing from {sdist.name}"
