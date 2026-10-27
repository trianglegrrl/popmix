"""The compiled module under test is built from the Rust sources as they are.

``maturin develop`` leaves the extension in ``python/popmix/``. One older
than the sources, or a second one under another ABI tag that Python
prefers (``_core.cpython-39-...so`` beside ``_core.abi3.so``), would test
old Rust code without saying so.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import popmix._core as core

ROOT = Path(__file__).parents[2]
PACKAGE = ROOT / "python" / "popmix"
EXTENSION = Path(core.__file__)

pytestmark = pytest.mark.skipif(PACKAGE not in EXTENSION.parents,
                                reason="popmix is installed from a wheel, not built here")


def test_one_compiled_module_in_the_package() -> None:
    built = sorted(p.name for p in PACKAGE.iterdir() if p.name.startswith("_core.")
                   and p.suffix in (".so", ".pyd"))
    assert built == [EXTENSION.name], "remove the others, then rebuild"


def test_the_compiled_module_is_newer_than_the_rust_sources() -> None:
    sources = [*(ROOT / "crates").rglob("*.rs"), *(ROOT / "crates").rglob("Cargo.toml"),
               ROOT / "Cargo.toml", ROOT / "Cargo.lock"]
    newest = max(sources, key=lambda p: p.stat().st_mtime)
    assert EXTENSION.stat().st_mtime >= newest.stat().st_mtime, (
        f"{newest.relative_to(ROOT)} changed after {EXTENSION.name} was built: rebuild")
