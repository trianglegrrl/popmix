#!/usr/bin/env python3
"""The AADR v66.p1 subset as a Zenodo archive, so the AADR goldens run anywhere.

`build_aadr_v66p1.py` builds the subset from the full release, which takes
the release itself and AdmixTools 8 convertf. This script packs that subset
(`aadr_v66p1.{geno,snp,ind}`, `manifest.json` and a README) into
`aadr_v66p1.tar.xz` for the Zenodo record, and fetches it back into the
reference cache, checking every file against the MD5 the build recorded:

    python tests/golden/aadr/archive_aadr_v66p1.py pack DIR
    python tests/golden/aadr/archive_aadr_v66p1.py fetch

AADR is released under CC0 1.0 (doi:10.7910/DVN/FFIDCW), and so is the
archive. `fetch` reads the record from $POPMIX_AADR_ZENODO_RECORD when set.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1"
NAME = "aadr_v66p1"
ARCHIVE = f"{NAME}.tar.xz"

# The subset's files as build_aadr_v66p1.py wrote them (its manifest.json).
MD5 = {"geno": "006c09f886d0039dd4beb3b8ba3cfe42",
       "snp": "50f66178fc81b8aa087cc4b135317e59",
       "ind": "2298265d7281b77fd65c7b7ece2e41ff"}

# The published Zenodo record and its archive; set when the record is published.
RECORD = os.environ.get("POPMIX_AADR_ZENODO_RECORD", "")
ARCHIVE_SHA256 = ""

README = """\
# AADR v66.p1 subset for the popmix golden tests

The populations the popmix AADR goldens use, taken from the Allen Ancient DNA
Resource v66.p1, 1240K public panel (Mallick et al. 2024,
doi:10.7910/DVN/FFIDCW), and relabelled `<label>.<data type>` (`French.DG`).
Sample IDs, SNPs (all 1,233,013, in order) and genotypes are v66.p1's own,
written as PACKEDANCESTRYMAP by AdmixTools 8 convertf. `manifest.json`
records the release's input checksums, the populations and the MD5 of each
file here.

v66.p1 is the patch release that withdrew the 161 genomes of Jacobs et al.
2019; this subset holds none of them.

Licence: CC0 1.0, as the AADR release
(https://creativecommons.org/publicdomain/zero/1.0/).

Built by tests/golden/aadr/build_aadr_v66p1.py and packed by
tests/golden/aadr/archive_aadr_v66p1.py in the popmix repository.
"""


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pack(src: Path, out: Path) -> Path:
    """Write `aadr_v66p1.tar.xz` from the built subset in `src` into `out`."""
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((src / "manifest.json").read_text())
    for entry in manifest.get("outputs", {}).values():
        entry["path"] = Path(entry["path"]).name
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
        (Path(tmp) / "README.md").write_text(README)
        target = out / ARCHIVE
        with tarfile.open(target, "w:xz") as tar:
            for ext in ("geno", "snp", "ind"):
                tar.add(src / f"{NAME}.{ext}", arcname=f"{NAME}.{ext}")
            tar.add(Path(tmp) / "manifest.json", arcname="manifest.json")
            tar.add(Path(tmp) / "README.md", arcname="README.md")
    return target


def installed(dest: Path, md5s: dict[str, str]) -> bool:
    """Whether every subset file in `dest` has its expected MD5."""
    return all((dest / f"{NAME}.{ext}").exists() and md5(dest / f"{NAME}.{ext}") == m for ext, m in md5s.items())


def install(archive: Path, dest: Path, md5s: dict[str, str]) -> None:
    """Unpack `archive` into `dest`, refusing unsafe members and wrong checksums.

    Raises
    ------
    ValueError
        If a member would land outside `dest`, is not a regular file, or a
        subset file's MD5 differs from `md5s`.
    """
    dest.mkdir(parents=True, exist_ok=True)
    root = dest.resolve()
    with tarfile.open(archive, "r:xz") as tar:
        members = tar.getmembers()
        for m in members:
            if not m.isfile() or root not in (root / m.name).resolve().parents:
                raise ValueError(f"unsafe archive member: {m.name}")
        with tempfile.TemporaryDirectory(dir=dest) as tmp:
            for m in members:
                with tar.extractfile(m) as f, open(Path(tmp) / m.name, "wb") as out:
                    shutil.copyfileobj(f, out)
            for ext, expected in md5s.items():
                got = md5(Path(tmp) / f"{NAME}.{ext}")
                if got != expected:
                    raise ValueError(f"{NAME}.{ext}: MD5 {got}, expected {expected}")
            for m in members:
                os.replace(Path(tmp) / m.name, dest / m.name)


def fetch(record: str = RECORD, dest: Path = CACHE) -> None:
    """Download the archive from Zenodo record `record` and install it in `dest`."""
    if installed(dest, MD5):
        print(f"{dest}: already installed")
        return
    if not record:
        raise SystemExit("No Zenodo record: set POPMIX_AADR_ZENODO_RECORD")
    url = f"https://zenodo.org/records/{record}/files/{ARCHIVE}?download=1"
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / ARCHIVE
        with urllib.request.urlopen(url) as response, open(target, "wb") as out:
            shutil.copyfileobj(response, out, 1 << 24)
        if ARCHIVE_SHA256 and sha256(target) != ARCHIVE_SHA256:
            raise ValueError(f"{url}: SHA-256 {sha256(target)}, expected {ARCHIVE_SHA256}")
        install(target, dest, MD5)
    print(f"{dest}: installed from {url}")


def main(argv: list[str]) -> None:
    if argv[:1] == ["pack"] and len(argv) == 2:
        if not installed(CACHE, MD5):
            raise SystemExit(f"{CACHE}: subset missing or changed; run build_aadr_v66p1.py")
        target = pack(CACHE, Path(argv[1]))
        print(f"{target}\nsha256 {sha256(target)}")
    elif argv == ["fetch"]:
        fetch()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
