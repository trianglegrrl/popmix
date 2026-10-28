#!/usr/bin/env python3
"""Record the environment each set of references was made in (TASKS.md K2).

Writes `environment.txt` into every tool-version directory under
tests/golden: the machine (`uname`), the tool's own version banner where it
prints one, the SHA-256 of each binary (several are patched, high-precision
builds named by their install path), and, for references made in R,
`sessionInfo()` with the packages loaded. Tools are run, never read.

Usage: python tests/golden/record_environment.py [directory...]
"""

from __future__ import annotations

import hashlib
import os
import platform
import subprocess
import sys
from pathlib import Path

GOLDEN = Path(__file__).resolve().parent
SHARE = Path.home() / ".local/share/popmix-references"
CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references"))
RSCRIPT = "/usr/bin/Rscript"
CONVERTF8 = SHARE / "admixtools-8.0.2/bin/convertf"

# Each tool-version directory: the R packages its generator loaded, the
# binaries it ran, a command printing the version, the Python package, or a
# file naming the commit.
TOOLS: dict[str, dict] = {
    "admixtools2/v2.0.10": {"r": ["admixtools"]},
    "formats/v2.0.10": {"r": ["admixtools"], "bins": [CONVERTF8]},
    "rlang/r-4.5.2": {"r": []},
    "admixture/v1.3.0": {"version": ["admixture", "--help"]},
    "plink2/v2.0.0-a.7LM": {"version": ["plink2", "--version"]},
    "dadi/v2.4.4": {"python": "dadi"},
    "dates/v753": {"bins": [CACHE / "dates-v753/bin/dates"]},
    "eigensoft/v7.2.1": {"bins": [SHARE / "eigensoft-7.2.1-hiprec/smartpca"]},
    "readv2/73739be": {"files": [SHARE / "readv2/COMMIT"]},
    "rolloff/v7.0.2": {"r": [], "bins": [SHARE / "admixtools-7.0.2/rolloff-hiprec"]},
}


def run(cmd: list[str]) -> str:
    out = subprocess.run(cmd, capture_output=True, text=True)
    return (out.stdout + out.stderr).strip()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binary(path: Path) -> str:
    resolved = Path(run(["which", str(path)]) or path)
    return f"{resolved}  sha256 {sha256(resolved)}"


def record(tool: dict) -> str:
    parts = [f"## machine\n{platform.platform()}\n{run(['uname', '-a'])}"]
    if "version" in tool:
        banner = [line for line in run(tool["version"]).splitlines() if "ersion" in line or "PLINK" in line]
        parts.append(f"## {' '.join(tool['version'])}\n{banner[0].strip('* ')}\n{binary(Path(tool['version'][0]))}")
    for path in tool.get("bins", []):
        parts.append(f"## binary\n{binary(Path(path))}")
    for path in tool.get("files", []):
        parts.append(f"## {path.name}\n{Path(path).read_text().strip()}")
    if "python" in tool:
        name = tool["python"]
        version = run([sys.executable, "-c", f"import importlib.metadata as m; print(m.version('{name}'))"])
        parts.append(f"## python package\n{name} {version}\nPython {platform.python_version()}")
    if "r" in tool:
        libs = "".join(f"suppressMessages(library({p}));" for p in tool["r"])
        parts.append(f"## R: sessionInfo()\n{run([RSCRIPT, '-e', libs + 'print(sessionInfo())'])}")
    return "\n\n".join(parts) + "\n"


def main(directories: list[str]) -> None:
    for directory in directories or sorted(TOOLS):
        (GOLDEN / directory / "environment.txt").write_text(record(TOOLS[directory]))
        print(GOLDEN / directory / "environment.txt")


if __name__ == "__main__":
    main(sys.argv[1:])
