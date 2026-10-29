"""Every licence the wheels must carry (TASKS.md C8): popmix's own, the BSD
notices of the code derived from other tools (NOTICE), and the licences of
the Rust crates compiled into the extension (THIRDPARTY.yml)."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
THIRDPARTY = yaml.safe_load((ROOT / "THIRDPARTY.yml").read_text())["third_party_libraries"]
NOTICE = (ROOT / "NOTICE").read_text()


def _license_files() -> str:
    return re.search(r"(?m)^license-files = \[(.*)\]", (ROOT / "pyproject.toml").read_text()).group(1)


@pytest.mark.parametrize("name", ["LICENSE", "NOTICE", "THIRDPARTY.yml"])
def test_the_wheels_carry_the_licence_file(name: str) -> None:
    assert f'"{name}"' in _license_files()


def test_every_crate_has_its_licence_text() -> None:
    assert [lib["package_name"] for lib in THIRDPARTY
            if not lib["licenses"] or any("NOT FOUND" in (x.get("text") or "NOT FOUND") for x in lib["licenses"])] == []


def test_every_crate_compiled_into_the_extension_is_listed() -> None:
    cargo = shutil.which("cargo")
    if cargo is None:
        pytest.skip("cargo not installed")
    meta = json.loads(subprocess.run([cargo, "metadata", "--format-version", "1", "--locked",
                                      "--manifest-path", str(ROOT / "crates/popmix-python/Cargo.toml")],
                                     capture_output=True, text=True, check=True).stdout)
    nodes = {n["id"]: n for n in meta["resolve"]["nodes"]}
    packages = {p["id"]: p for p in meta["packages"]}
    root = next(p["id"] for p in meta["packages"] if p["name"] == "popmix-python")
    seen, todo = set(), [root]
    while todo:
        for dep in nodes[todo.pop()]["deps"]:
            normal = any(k["kind"] is None for k in dep["dep_kinds"])
            if normal and dep["pkg"] not in seen:
                seen.add(dep["pkg"])
                todo.append(dep["pkg"])
    compiled = {(packages[i]["name"], packages[i]["version"]) for i in seen if packages[i]["source"]}
    listed = {(lib["package_name"], lib["package_version"]) for lib in THIRDPARTY}
    assert sorted(compiled - listed) == [], "run .github/scripts/thirdparty-licenses.py"


@pytest.mark.parametrize("clause", [
    # each component's own licence, verbatim
    "Neither the name Broad Institute, Inc. nor the names of its",  # ADMIXTOOLS
    "Neither the name Broad Institute, Inc. Harvard University, nor the names of its",  # EIGENSOFT
    "Neither the name of the <organization> nor the",  # L-BFGS-B 3.0
    "3. The names of its contributors may not be used to endorse or promote",  # mt19937ar.c
])
def test_notice_quotes_each_licence(clause: str) -> None:
    assert clause in NOTICE


def test_notice_names_every_l_bfgs_b_author() -> None:
    authors = ["Ciyou Zhu", "Richard Byrd", "Peihuang Lu", "Jorge Nocedal", "Jose Luis Morales"]
    assert [a for a in authors if a not in NOTICE] == []


@pytest.mark.parametrize("name", ["LICENSE", "NOTICE", "THIRDPARTY.yml"])
def test_the_release_checks_each_built_wheel_carries_the_licence_file(name: str) -> None:
    check = (ROOT / ".github/scripts/check-wheel.sh").read_text()
    assert f"licenses/{name}" in check
