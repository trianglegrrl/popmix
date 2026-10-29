"""The commands take and give the same things the same way (TASKS.md M4).

A table goes to stdout (or ``-o``) alone, with notes on stderr; scans take
``--format tsv|json``; an option of the Python function is an option of
its command; JSON ends with a newline; "all populations" leaves out
``Ignore``; the docs list every format convertf reads.
"""

from __future__ import annotations

import io
import json
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")
DG6 = str(ROOT / "tests" / "golden" / "admixtools2" / "v2.0.10" / "dg_6pop" / "f2_dir")
QPADM = ["--target", "French.DG", "--left", "Han.DG,Papuan.DG", "--right", "Mbuti.DG,Yoruba.DG,Chimp.REF"]


def popmix(*args: str) -> subprocess.CompletedProcess:
    out = subprocess.run(["popmix", *args], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out


# --- sfs -----------------------------------------------------------------------

def test_sfs_of_one_population_writes_only_the_table_to_stdout():
    out = popmix("sfs", "-i", SMALL, "--populations", "French")
    table = pd.read_csv(io.StringIO(out.stdout), sep="\t")
    assert list(table.columns) == ["French", "n", "masked"]
    assert "tajima_d" in out.stderr


def test_sfs_summary_stays_off_the_table_file(tmp_path):
    path = tmp_path / "sfs.tsv"
    out = popmix("sfs", "-i", SMALL, "--populations", "French", "-o", str(path))
    assert out.stdout == ""
    assert path.read_text(encoding="utf-8").startswith("French\tn\tmasked\n")


@pytest.mark.parametrize("rule", ["detect", "admixtools2", "diploid"])
def test_sfs_takes_pseudohaploid(rule):
    expected = pm.sfs(SMALL, ["French"], pseudohaploid=rule).data.ravel()
    out = popmix("sfs", "-i", SMALL, "--populations", "French", "--pseudohaploid", rule).stdout
    assert list(pd.read_csv(io.StringIO(out), sep="\t")["n"]) == pytest.approx(list(expected), rel=1e-12)


@pytest.mark.parametrize("args", [
    ["sfs", "-i", SMALL, "--populations", "French"],
    ["kinship", "-i", SMALL, "--populations", "French,Yoruba"],
    ["kinship", "-i", SMALL, "--populations", "French,Yoruba", "--method", "read"],
], ids=["sfs", "king", "read"])
def test_scans_write_json_with_format_json(args):
    tsv = pd.read_csv(io.StringIO(popmix(*args).stdout), sep="\t")
    text = popmix(*args, "--format", "json").stdout
    assert text.endswith("\n")
    rows = json.loads(text)
    assert len(rows) == len(tsv)
    assert list(rows[0]) == list(tsv.columns)


@pytest.fixture(scope="module")
def dg6():
    return pm.read_admixtools2_f2(DG6)


@pytest.fixture(scope="module")
def f2_dir(dg6, tmp_path_factory) -> str:
    path = str(tmp_path_factory.mktemp("dg6") / "f2")
    dg6.save(path)
    return path


def test_json_tables_end_with_a_newline(f2_dir):
    out = popmix("qpdstat", "--f2-dir", f2_dir, "--pop1", "Mbuti.DG", "--pop2", "Han.DG", "--pop3", "French.DG",
                 "--pop4", "Papuan.DG", "--format", "json")
    assert out.stdout.endswith("\n")
    assert len(json.loads(out.stdout)) == 1


# --- qpadm and qpwave ----------------------------------------------------------

def test_qpadm_prints_weight_ses_in_text_and_tsv(dg6, f2_dir):
    r = pm.qpadm(dg6, ["Han.DG", "Papuan.DG"], ["Mbuti.DG", "Yoruba.DG", "Chimp.REF"], "French.DG")
    text = popmix("qpadm", "--f2-dir", f2_dir, *QPADM).stdout
    for se in r.weight_se:
        assert f"{se:.4f}" in text
    tsv = pd.read_csv(io.StringIO(popmix("qpadm", "--f2-dir", f2_dir, *QPADM, "--format", "tsv").stdout), sep="\t")
    assert list(tsv["se"]) == pytest.approx(list(r.weight_se), abs=1e-4)


def test_qpadm_takes_fudge(dg6, f2_dir):
    r = pm.qpadm(dg6, ["Han.DG", "Papuan.DG"], ["Mbuti.DG", "Yoruba.DG", "Chimp.REF"], "French.DG", fudge=0.1)
    out = json.loads(popmix("qpadm", "--f2-dir", f2_dir, *QPADM, "--fudge", "0.1", "--format", "json").stdout)
    assert out["chisq"] == pytest.approx(r.chisq, rel=1e-12)
    default = json.loads(popmix("qpadm", "--f2-dir", f2_dir, *QPADM, "--format", "json").stdout)
    assert default["chisq"] != pytest.approx(r.chisq, rel=1e-6)


def test_qpwave_takes_fudge(dg6, f2_dir):
    r = pm.qpwave(dg6, ["French.DG", "Han.DG", "Papuan.DG"], ["Mbuti.DG", "Yoruba.DG", "Chimp.REF"], fudge=0.1)
    out = json.loads(popmix("qpwave", "--f2-dir", f2_dir, "--left", "French.DG,Han.DG,Papuan.DG",
                            "--right", "Mbuti.DG,Yoruba.DG,Chimp.REF", "--fudge", "0.1", "--format", "json").stdout)
    assert [t["chisq"] for t in out["rank_tests"]] == pytest.approx([t.chisq for t in r.rank_tests], rel=1e-12)


# --- dates and rolloff ---------------------------------------------------------

DATES = ["dates", "-i", SMALL, "--admixed", "Simulation.1", "--sources", "French,Yoruba"]


def test_dates_takes_generation_time_and_no_jackknife():
    r = pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], generation_time=25.0, jackknife=False)
    out = json.loads(popmix(*DATES, "--generation-time", "25", "--no-jackknife", "--format", "json").stdout)
    assert out["generations"] == pytest.approx(r.generations, rel=1e-12)
    assert out["years_bp"] == pytest.approx(r.years_bp, rel=1e-12)
    assert out["years_bp"] == pytest.approx(25 * r.generations, rel=1e-12)
    text = popmix(*DATES, "--generation-time", "25", "--no-jackknife").stdout
    assert "25 years/gen" in text


def test_dates_json_names_its_populations_as_its_options():
    out = json.loads(popmix(*DATES, "--no-jackknife", "--format", "json").stdout)
    assert out["admixed"] == "Simulation.1"
    assert out["sources"] == ["French", "Yoruba"]


def test_rolloff_takes_its_fit_range():
    r = pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"], lval=1.0, hval=20.0)
    out = popmix("rolloff", "-i", SMALL, "--admixed", "Simulation.1", "--sources", "French,Yoruba",
                 "--lval", "1", "--hval", "20").stdout
    assert f"decay {r.decay:.3f} generations" in out


# --- Ignore --------------------------------------------------------------------

def test_extract_f2_of_all_populations_leaves_out_ignore(tmp_path):
    prefix = tmp_path / "data"
    for ext in ("geno", "snp"):
        shutil.copy(f"{SMALL}.{ext}", f"{prefix}.{ext}")
    lines = Path(f"{SMALL}.ind").read_text(encoding="utf-8").splitlines()
    lines[0] = lines[0].rsplit(" ", 1)[0] + " Ignore"
    Path(f"{prefix}.ind").write_text("\n".join(lines) + "\n", encoding="utf-8")
    popmix("extract-f2", "-i", str(prefix), "-o", str(tmp_path / "f2"))
    assert "Ignore" not in pm.read_f2(str(tmp_path / "f2")).populations


# --- formats in the docs -------------------------------------------------------

def _input_formats() -> list[str]:
    from popmix.cli import cli

    option = next(p for p in cli.commands["convertf"].params if p.name == "input_format")
    return list(option.type.choices)


def test_convertf_help_names_every_input_format():
    text = popmix("convertf", "--help").stdout.lower()
    assert "6 input formats" not in text
    listed = {name for item in re.findall(r"^\s*- ([a-z/]+):", text, re.MULTILINE) for name in item.split("/")}
    assert listed == set(_input_formats())


def test_readme_names_every_input_format():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "reads 6 " not in readme
    inputs = next(line for line in readme.splitlines() if line.startswith("**Input formats:**"))
    for fmt in _input_formats():
        assert f"`{fmt}`" in inputs, fmt
    assert "### TGENO" in readme


# --- compat --------------------------------------------------------------------

# the command of each public function that takes compat
COMPAT_COMMANDS = {"qpadm": "qpadm", "dates": "dates", "smartpca": "pca", "read_kinship": "kinship",
                   "joint_sfs": "sfs", "write_dot": "qpgraph", "qpgraph": "qpgraph"}


def test_every_compat_switch_is_on_the_cli():
    """DEVIATIONS.md step 2: Python and the CLI both take the compat switch."""
    import inspect

    from popmix.cli import cli

    with_compat = {name for name in pm.__all__ if callable(getattr(pm, name))
                   and not isinstance(getattr(pm, name), type)
                   and "compat" in inspect.signature(getattr(pm, name)).parameters}
    assert with_compat == set(COMPAT_COMMANDS)
    for command in COMPAT_COMMANDS.values():
        assert "compat" in {p.name for p in cli.commands[command].params}, command


@pytest.mark.parametrize("compat", [None, "admixtools2-2.0.10"])
def test_qpgraph_dot_takes_compat(dg6, f2_dir, tmp_path, compat):
    """--dot writes what write_dot does, admixtools2's text with --compat (ADMIXTOOLS2-0004)."""
    ref = json.loads((ROOT / "tests/golden/admixtools2/v2.0.10/dg_6pop/qpgraph.json").read_text(encoding="utf-8"))
    edges = next(g for g in ref["graphs"] if g["name"] == "admix1")["edges"]
    graph = tmp_path / "graph.txt"
    graph.write_text("".join(f"{a}\t{b}\n" for a, b in edges), encoding="utf-8")
    dot = tmp_path / "graph.dot"
    popmix("qpgraph", "--f2-dir", f2_dir, "--graph", str(graph), "--dot", str(dot),
           *(["--compat", compat] if compat else []))
    expected = pm.write_dot(pm.qpgraph(dg6, [tuple(e) for e in edges]), compat=compat)
    assert dot.read_text(encoding="utf-8") == expected


def test_qpgraph_compat_fits_as_admixtools2(dg6, f2_dir, tmp_path):
    """--compat stops the fit where admixtools2 does (ADMIXTOOLS2-0009), with or without --dot."""
    ref = json.loads((ROOT / "tests/golden/admixtools2/v2.0.10/dg_6pop/qpgraph.json").read_text(encoding="utf-8"))
    edges = next(g for g in ref["graphs"] if g["name"] == "admix1")["edges"]
    graph = tmp_path / "graph.txt"
    graph.write_text("".join(f"{a}\t{b}\n" for a, b in edges), encoding="utf-8")
    out = json.loads(popmix("qpgraph", "--f2-dir", f2_dir, "--graph", str(graph), "--compat", "admixtools2-2.0.10",
                            "--format", "json").stdout)
    expected = pm.qpgraph(dg6, [tuple(e) for e in edges], compat="admixtools2-2.0.10")
    assert out["score"] == pytest.approx(expected.score, rel=1e-12)
