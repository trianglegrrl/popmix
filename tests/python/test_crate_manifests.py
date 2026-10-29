"""The Rust crates can be published, and builds are reproducible (TASKS.md L1).

popmix-core and popmix-io are libraries others may depend on from
crates.io: their path dependencies carry the version too, and every crate
states the oldest Rust it supports. popmix-python is only the extension
module inside the wheel and is never published on its own. Cargo.lock is
tracked, so the wheel is built from the dependency versions that were
tested.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    tomllib = pytest.importorskip("tomli")

ROOT = Path(__file__).parents[2]
CRATES = sorted(p.parent for p in (ROOT / "crates").glob("*/Cargo.toml"))


def _toml(path: Path) -> dict:
    return tomllib.loads(path.read_text(encoding="utf-8"))


WORKSPACE = _toml(ROOT / "Cargo.toml")["workspace"]


@pytest.mark.parametrize("crate", CRATES, ids=lambda p: p.name)
def test_path_dependencies_carry_the_workspace_version(crate: Path) -> None:
    deps = _toml(crate / "Cargo.toml").get("dependencies", {})
    for name, spec in deps.items():
        if isinstance(spec, dict) and "path" in spec:
            assert spec.get("version") == f"={WORKSPACE['package']['version']}", name


@pytest.mark.parametrize("crate", CRATES, ids=lambda p: p.name)
def test_every_crate_states_the_rust_it_needs(crate: Path) -> None:
    channel = _toml(ROOT / "rust-toolchain.toml")["toolchain"]["channel"]
    assert WORKSPACE["package"]["rust-version"] == channel
    assert _toml(crate / "Cargo.toml")["package"]["rust-version"] == {"workspace": True}


def test_only_the_libraries_are_publishable() -> None:
    publish = {c.name: _toml(c / "Cargo.toml")["package"].get("publish", True) for c in CRATES}
    assert publish == {"popmix-core": True, "popmix-io": True, "popmix-python": False}


def test_the_lockfile_is_tracked() -> None:
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", "Cargo.lock"], cwd=ROOT,
                             capture_output=True, text=True)
    assert tracked.returncode == 0, tracked.stderr


@pytest.mark.skipif(shutil.which("cargo") is None, reason="cargo not available")
@pytest.mark.parametrize("crate", ["popmix-core", "popmix-io"])
def test_the_published_crate_carries_the_licence_and_notices(crate: str) -> None:
    listed = subprocess.run(["cargo", "package", "-p", crate, "--list", "--allow-dirty"], cwd=ROOT,
                            capture_output=True, text=True, check=True).stdout.split()
    for name in ("LICENSE", "NOTICE"):
        assert name in listed, name
    # what the packaged files say is the repository's own text
    for name in ("LICENSE", "NOTICE"):
        assert (ROOT / "crates" / crate / name).read_text(encoding="utf-8") == (ROOT / name).read_text(encoding="utf-8")
