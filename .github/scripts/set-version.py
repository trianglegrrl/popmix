#!/usr/bin/env python3
"""Set the version a release build carries from the ref it runs on.

    python .github/scripts/set-version.py REF [ROOT]

A tag `v1.0.0` must name the version in Cargo.toml and pyproject.toml, and
changes nothing. A test tag `v1.0.0rc1` builds 1.0.0rc1: pyproject.toml gets
`1.0.0rc1`, and the workspace version and the pinned path dependencies get
the Cargo form `1.0.0-rc.1`; release.yml sends those builds to TestPyPI
only. Any ref not starting with `v` (a branch, run by hand) changes nothing.
Prints the Python version the build carries.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

CRATES = ["crates/popmix-python/Cargo.toml", "crates/popmix-io/Cargo.toml"]


def main(ref: str, root: Path) -> None:
    workspace = root / "Cargo.toml"
    version = re.search(r'(?m)^version = "([^"]+)"', workspace.read_text()).group(1)
    if not ref.startswith("v"):
        print(version)
        return
    tag = re.fullmatch(r"v(\d+\.\d+\.\d+)(?:rc(\d+))?", ref)
    if not tag:
        sys.exit(f"{ref}: a release tag is v<version> or v<version>rc<n>")
    if tag.group(1) != version:
        sys.exit(f"{ref}: the tag names {tag.group(1)}, but Cargo.toml and pyproject.toml hold {version}")
    if tag.group(2) is None:
        print(version)
        return
    python, cargo = f"{version}rc{tag.group(2)}", f"{version}-rc.{tag.group(2)}"
    for path, old, new in [(workspace, f'version = "{version}"', f'version = "{cargo}"'),
                           (root / "pyproject.toml", f'version = "{version}"', f'version = "{python}"')]:
        text = path.read_text()
        path.write_text(text.replace(old, new, 1))
    for name in CRATES:
        path = root / name
        path.write_text(path.read_text().replace(f'version = "={version}"', f'version = "={cargo}"'))
    print(python)


if __name__ == "__main__":
    main(sys.argv[1], Path(sys.argv[2]) if len(sys.argv) > 2 else Path("."))
