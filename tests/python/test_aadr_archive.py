"""The AADR v66.p1 subset archive: packed, fetched and verified (TASKS.md K2).

Runs on small stand-in files, with no network: `pack` writes the archive the
Zenodo record holds, `install` unpacks it and checks every file against its
expected MD5.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import tarfile
from io import BytesIO
from pathlib import Path

import pytest

AADR_DIR = Path(__file__).resolve().parents[1] / "golden" / "aadr"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, AADR_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


archive = load("archive_aadr_v66p1")


def subset(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    """A stand-in for the built subset: three files and a manifest."""
    src = tmp_path / "cache"
    src.mkdir()
    md5s = {}
    for ext, body in {"geno": b"GENO  data", "snp": b"rs1 1 0.0 100 A G\n", "ind": b"S1 M French.DG\n"}.items():
        (src / f"aadr_v66p1.{ext}").write_bytes(body)
        md5s[ext] = hashlib.md5(body).hexdigest()
    manifest = {"release": "AADR v66.p1", "outputs": {ext: {"path": f"/home/x/aadr_v66p1.{ext}", "md5": m}
                                                      for ext, m in md5s.items()}}
    (src / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return src, md5s


def test_pack_then_install_gives_the_same_files(tmp_path):
    src, md5s = subset(tmp_path)
    packed = archive.pack(src, tmp_path / "out")
    dest = tmp_path / "dest"
    archive.install(packed, dest, md5s)
    assert (dest / "aadr_v66p1.geno").read_bytes() == (src / "aadr_v66p1.geno").read_bytes()


def test_the_packed_manifest_names_files_not_local_paths(tmp_path):
    src, md5s = subset(tmp_path)
    dest = tmp_path / "dest"
    archive.install(archive.pack(src, tmp_path / "out"), dest, md5s)
    manifest = json.loads((dest / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["outputs"]["geno"]["path"] == "aadr_v66p1.geno"


def test_the_archive_carries_a_readme_with_the_licence(tmp_path):
    src, md5s = subset(tmp_path)
    dest = tmp_path / "dest"
    archive.install(archive.pack(src, tmp_path / "out"), dest, md5s)
    assert "CC0 1.0" in (dest / "README.md").read_text(encoding="utf-8")


def test_install_refuses_a_file_with_the_wrong_checksum(tmp_path):
    src, md5s = subset(tmp_path)
    packed = archive.pack(src, tmp_path / "out")
    with pytest.raises(ValueError, match="aadr_v66p1.snp"):
        archive.install(packed, tmp_path / "dest", {**md5s, "snp": "0" * 32})


def test_install_refuses_a_member_outside_the_destination(tmp_path):
    bad = tmp_path / "bad.tar.xz"
    with tarfile.open(bad, "w:xz") as tar:
        info = tarfile.TarInfo("../escaped")
        info.size = 1
        tar.addfile(info, BytesIO(b"x"))
    with pytest.raises(ValueError, match="escaped"):
        archive.install(bad, tmp_path / "dest", {})
    assert not (tmp_path / "escaped").exists()


def test_installed_files_are_kept_when_they_already_match(tmp_path):
    src, md5s = subset(tmp_path)
    dest = tmp_path / "dest"
    archive.install(archive.pack(src, tmp_path / "out"), dest, md5s)
    assert archive.installed(dest, md5s) is True


def test_a_changed_installed_file_does_not_count_as_installed(tmp_path):
    src, md5s = subset(tmp_path)
    dest = tmp_path / "dest"
    archive.install(archive.pack(src, tmp_path / "out"), dest, md5s)
    (dest / "aadr_v66p1.ind").write_bytes(b"changed")
    assert archive.installed(dest, md5s) is False


def test_the_expected_checksums_are_those_the_build_recorded():
    manifest = archive.CACHE / "manifest.json"
    if not manifest.exists():
        pytest.skip(f"AADR v66.p1 subset not built: {manifest}")
    outputs = json.loads(manifest.read_text(encoding="utf-8"))["outputs"]
    assert {ext: entry["md5"] for ext, entry in outputs.items()} == archive.MD5


def test_fetch_reads_the_published_record_by_default(monkeypatch):
    monkeypatch.delenv("POPMIX_AADR_ZENODO_RECORD", raising=False)
    assert load("archive_aadr_v66p1").RECORD == "23750150"


def test_fetch_refuses_an_archive_with_the_wrong_sha256(tmp_path, monkeypatch):
    monkeypatch.setattr(archive.urllib.request, "urlopen", lambda url: BytesIO(b"not the archive"))
    with pytest.raises(ValueError, match="SHA-256"):
        archive.fetch("23750150", tmp_path / "dest")
