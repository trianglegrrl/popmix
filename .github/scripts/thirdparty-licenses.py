#!/usr/bin/env python3
"""Write THIRDPARTY.yml: the licence of every Rust crate compiled into the
extension, which the wheels carry (TASKS.md C8).

    python .github/scripts/thirdparty-licenses.py

cargo-bundle-licenses (`cargo install --locked cargo-bundle-licenses`) reads
each crate's licence files. A few crates publish none; for those, MIT (each
is MIT or offers it) is filled in from the standard text, with the crate's
authors as its holders. Rerun after changing Cargo.lock;
tests/python/test_thirdparty_licenses.py fails until then.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "crates/popmix-python/Cargo.toml"
OUT = ROOT / "THIRDPARTY.yml"

MIT = """\
MIT License

Copyright (c) {holders}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


def missing(entry: dict) -> bool:
    return "NOT FOUND" in (entry.get("text") or "NOT FOUND")


def main() -> None:
    raw = subprocess.run(["cargo-bundle-licenses", "--format", "yaml", "--output", "-"],
                         cwd=MANIFEST.parent, capture_output=True, text=True, check=True).stdout
    bundle = yaml.safe_load(raw)
    meta = json.loads(subprocess.run(["cargo", "metadata", "--format-version", "1", "--locked",
                                      "--manifest-path", str(MANIFEST)],
                                     capture_output=True, text=True, check=True).stdout)
    authors = {(p["name"], p["version"]): p["authors"] for p in meta["packages"]}
    for lib in bundle["third_party_libraries"]:
        found = [x for x in lib["licenses"] if not missing(x)]
        if found:
            lib["licenses"] = found
            continue
        if "MIT" not in lib["license"].replace("(", " ").replace(")", " ").split():
            sys.exit(f"{lib['package_name']} {lib['package_version']}: no licence text, and not MIT ({lib['license']})")
        names = authors[(lib["package_name"], lib["package_version"])]
        holders = ", ".join(n.split(" <")[0] for n in names) or f"the {lib['package_name']} authors"
        lib["licenses"] = [{"license": "MIT", "text": MIT.format(holders=holders)}]
    OUT.write_text(yaml.safe_dump(bundle, sort_keys=False, allow_unicode=True, width=1000), encoding="utf-8")
    print(f"{OUT}: {len(bundle['third_party_libraries'])} crates")


if __name__ == "__main__":
    main()
