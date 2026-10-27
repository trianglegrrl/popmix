"""Every Rust dependency a crate declares is used by that crate (TASKS.md L2).

A dependency counts as used when its library name (``-`` read as ``_``)
appears as a path (``name::``) or an extern crate in the crate's sources,
including its tests and benches.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    tomllib = pytest.importorskip("tomli")

ROOT = Path(__file__).parents[2]
CRATES = sorted(p.parent for p in (ROOT / "crates").glob("*/Cargo.toml"))


def _declared(manifest: dict) -> list[str]:
    return [name for table in ("dependencies", "dev-dependencies")
            for name in manifest.get(table, {})]


@pytest.mark.parametrize("crate", CRATES, ids=lambda p: p.name)
def test_every_declared_dependency_is_used(crate: Path) -> None:
    manifest = tomllib.loads((crate / "Cargo.toml").read_text())
    source = "\n".join(p.read_text() for p in crate.rglob("*.rs") if "target" not in p.parts)
    unused = [name for name in _declared(manifest)
              if not re.search(rf"\b{name.replace('-', '_')}(::|\s*;|\s+as\b)", source)]
    assert unused == []


def test_every_workspace_dependency_is_used_by_a_crate() -> None:
    workspace = tomllib.loads((ROOT / "Cargo.toml").read_text())["workspace"]["dependencies"]
    declared = {name for crate in CRATES
                for name in _declared(tomllib.loads((crate / "Cargo.toml").read_text()))}
    assert sorted(set(workspace) - declared) == []


def test_public_rust_types_are_ndarray_not_nalgebra_or_faer() -> None:
    """TASKS.md I5: ndarray is the array type of the public API; nalgebra
    (small dense solves) and faer (PCA's eigendecomposition) stay inside."""
    internal = re.compile(r"\b(DMatrix|DVector|nalgebra::|faer::)")
    public = re.compile(r"^\s*pub (fn|struct|enum|type|const|static|[a-z_][a-z0-9_]*\s*:)")
    found = []
    for crate in ("popmix-core", "popmix-io"):
        for path in sorted((ROOT / "crates" / crate / "src").rglob("*.rs")):
            if path.name == "tests.rs" or "tests" in path.parts:
                continue
            signature = None
            for number, line in enumerate(path.read_text().splitlines(), 1):
                if public.match(line):
                    signature = [line]
                elif signature is not None and not line.rstrip().endswith(("{", ";", ",")):
                    signature.append(line)
                if signature is not None and internal.search(" ".join(signature)):
                    found.append(f"{path.relative_to(ROOT)}:{number}")
                if signature is not None and line.rstrip().endswith(("{", ";", ",")):
                    signature = None
    assert found == []
