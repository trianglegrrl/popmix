"""The version and project URLs agree everywhere a release reads them.

PyPI reads pyproject.toml, the compiled module reads the Cargo workspace
version, and bioconda reads recipes/meta.yaml.
"""

import re
import shutil
import subprocess
import sys
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



LICENSE = "MIT AND BSD-3-Clause"


def test_the_licence_agrees_everywhere() -> None:
    """popmix is MIT; the code derived from ADMIXTOOLS and EIGENSOFT keeps their
    BSD-3 licence, whose notices are in NOTICE."""
    pyproject = _read("pyproject.toml")
    assert _match(pyproject, r'^license = "([^"]+)"') == LICENSE
    assert "NOTICE" in _match(pyproject, r"^license-files = \[(.*)\]")
    assert _match(_read("Cargo.toml"), r'^license = "([^"]+)"') == LICENSE
    recipe = _read("recipes/meta.yaml")
    assert _match(recipe, r"^  license: (.+)$") == LICENSE
    assert "    - NOTICE" in recipe
    notice = _read("NOTICE")
    for holder in ("Broad Institute, Inc.", "Harvard Medical School"):
        assert holder in notice


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
        "NOTICE",
        "README.md",
        "Cargo.lock",
        "pyproject.toml",
        "crates/popmix-core/Cargo.toml",
        "crates/popmix-io/Cargo.toml",
        "crates/popmix-python/Cargo.toml",
        "python/popmix/__init__.py",
    ):
        assert needed in names, f"{needed} missing from {sdist.name}"


def _requirement_names(block: str) -> set[str]:
    names = re.findall(r'^\s*"([A-Za-z0-9_.-]+)', block, re.MULTILINE)
    return {n.lower() for n in names if n.lower() != "popmix"}


def test_cli_starts_with_only_the_required_dependencies() -> None:
    """`pip install popmix` must give a working `popmix` command: block every
    package that only an optional extra provides, then start the CLI."""
    text = _read("pyproject.toml")
    core = _requirement_names(text.split("\ndependencies = [")[1].split("]")[0])
    extras = text.split("[project.optional-dependencies]")[1].split("\n[project.")[0]
    groups = re.split(r"^\w+ = \[", extras, flags=re.MULTILINE)
    optional = set().union(*map(_requirement_names, groups)) - core
    # import names that differ from the distribution names
    blocked = sorted({"pyyaml": "yaml"}.get(n, n).replace("-", "_") for n in optional)
    script = f"""
import importlib.abc, sys
BLOCKED = {blocked!r}
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name.split(".")[0] in BLOCKED:
            raise ModuleNotFoundError(f"blocked optional dependency {{name}}")
sys.meta_path.insert(0, Block())
from popmix.cli import main
for args in (["--version"], ["--help"]):
    sys.argv = ["popmix", *args]
    try:
        main()
    except SystemExit as e:
        assert e.code == 0, (args, e.code)
"""
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_only_the_python_package_and_cli_promise_a_stable_api() -> None:
    """1.0 is the Python package and the CLI; the Rust crates are its
    internals, without a stability promise (TASKS.md I6)."""
    readme = _read("README.md")
    stability = " ".join(readme.split("## API stability", 1)[1].split("\n## ", 1)[0].split())
    assert "Rust crates" in stability
    assert "no stability promise" in stability
    for crate in ("popmix-core", "popmix-io"):
        doc = " ".join(line[3:].strip() for line in _read(f"crates/{crate}/src/lib.rs").splitlines()
                       if line.startswith("//!"))
        assert "no stability promise" in doc, crate
