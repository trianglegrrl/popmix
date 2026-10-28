"""A release tag sets the version the wheels carry (TASKS.md K8).

`v1.0.0` must match the version in Cargo.toml and pyproject.toml; a test tag
`v1.0.0rc1` builds 1.0.0rc1 (Cargo `1.0.0-rc.1`), which goes to TestPyPI
only. Any other ref (a branch, run by hand) builds the version as it is.
"""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / ".github/scripts/set-version.py"
FILES = ["Cargo.toml", "pyproject.toml", "crates/popmix-python/Cargo.toml", "crates/popmix-io/Cargo.toml"]


def copy(tmp_path: Path) -> Path:
    for name in FILES:
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / name, tmp_path / name)
    return tmp_path


def run(ref: str, root: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), ref, str(root)], capture_output=True, text=True)


def test_a_final_tag_matching_the_version_changes_nothing(tmp_path):
    root = copy(tmp_path)
    assert run("v1.0.0", root).returncode == 0
    assert [(root / f).read_text() for f in FILES] == [(ROOT / f).read_text() for f in FILES]


def test_a_branch_changes_nothing(tmp_path):
    root = copy(tmp_path)
    assert run("main", root).returncode == 0
    assert (root / "pyproject.toml").read_text() == (ROOT / "pyproject.toml").read_text()


def test_a_test_tag_sets_the_python_version(tmp_path):
    root = copy(tmp_path)
    run("v1.0.0rc1", root)
    assert 'version = "1.0.0rc1"' in (root / "pyproject.toml").read_text()


def test_a_test_tag_sets_the_workspace_version(tmp_path):
    root = copy(tmp_path)
    run("v1.0.0rc1", root)
    assert 'version = "1.0.0-rc.1"' in (root / "Cargo.toml").read_text()


def test_a_test_tag_moves_the_pinned_path_dependencies(tmp_path):
    root = copy(tmp_path)
    run("v1.0.0rc2", root)
    for name in ("crates/popmix-python/Cargo.toml", "crates/popmix-io/Cargo.toml"):
        text = (root / name).read_text()
        assert 'version = "=1.0.0-rc.2"' in text, name
        assert 'version = "=1.0.0"' not in text, name


def test_a_tag_for_another_version_is_refused(tmp_path):
    out = run("v1.0.1", copy(tmp_path))
    assert out.returncode != 0
    assert "1.0.0" in out.stderr


def test_a_tag_of_another_shape_is_refused(tmp_path):
    assert run("v1.0.0-beta", copy(tmp_path)).returncode != 0


def test_prints_the_python_version_it_builds(tmp_path):
    assert run("v1.0.0rc3", copy(tmp_path)).stdout.strip() == "1.0.0rc3"
