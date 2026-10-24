"""Run the popmix nf-core modules in Nextflow and check what they produce.

The modules live in an nf-core/modules checkout (NFCORE_MODULES, default
~/nf-core/modules). A small workflow including them is written and run with
the local executor and the popmix on PATH:

- on AADR, extract-f2 for the dg_6pop populations, then f2, f3, f4, qpwave,
  qpadm and qpgraph on its f2 blocks, compared with admixtools2 2.0.10 at the
  golden-test tolerances;
- on the small DATES test data, findgraphs, dates, convertf and merge,
  compared with popmix run directly.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

import popmix as pm

REPO = Path(__file__).resolve().parents[2]
MODULES = Path(os.environ.get("NFCORE_MODULES", "~/nf-core/modules")).expanduser() / "modules" / "nf-core" / "popmix"
AADR = "/references/aadr/v62.0_1240k_public"
DG6 = REPO / "tests" / "golden" / "admixtools2" / "v2.0.10" / "dg_6pop"
SMALL = REPO / "tests" / "golden" / "dates" / "family_chr20_22"
FIXTURES = REPO / "tests" / "fixtures"
SMALL_POPS = ["Yoruba", "French", "Simulation", "Simulation.1"]

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not installed"),
    pytest.mark.skipif(not MODULES.is_dir(), reason=f"no popmix modules in {MODULES}"),
    pytest.mark.skipif(not os.path.exists(f"{AADR}.geno"), reason="AADR data not available"),
]


def _ref(file: str) -> dict:
    return json.loads((DG6 / file).read_text())


def _groovy(value) -> str:
    if isinstance(value, list):
        return "[" + ", ".join(_groovy(v) for v in value) + "]"
    return "'" + str(value).replace("\\", "\\\\").replace("'", "\\'") + "'"


def _files(prefix) -> str:
    return ", ".join(f"file('{prefix}.{ext}', checkIfExists: true)" for ext in ("geno", "snp", "ind"))


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """Run the workflow once; return a function giving a module's published outputs."""
    work = tmp_path_factory.mktemp("nextflow")
    f2_row = _ref("f2.json")["f2"][0]
    f3_row = _ref("f3_f4.json")["f3"][0]
    f4_row = _ref("f3_f4.json")["f4"][0]
    qpwave = _ref("qpwave.json")["models"][0]
    qpadm = next(m for m in _ref("qpadm.json")["models"] if not m["constrained"])
    graph = next(g for g in _ref("qpgraph.json")["graphs"] if g["name"] == "admix1")
    graph_file = work / "admix1.txt"
    graph_file.write_text("".join(f"{a}\t{b}\n" for a, b in graph["edges"]))

    names = ["extractf2", "f2", "f3", "f4", "qpwave", "qpadm", "qpgraph", "findgraphs", "dates", "convertf",
             "merge"]
    includes = "\n".join(f"include {{ POPMIX_{n.upper()} }} from '{MODULES / n / 'main'}'" for n in names)
    includes += f"\ninclude {{ POPMIX_EXTRACTF2 as POPMIX_EXTRACTF2_SMALL }} from '{MODULES / 'extractf2' / 'main'}'"
    (work / "main.nf").write_text(f"""
{includes}

workflow {{
    POPMIX_EXTRACTF2([ [ id:'dg_6pop' ], {_files(AADR)} ], {_groovy(_ref("f2.json")["pops"])})
    f2 = POPMIX_EXTRACTF2.out.f2_blocks
    POPMIX_F2(f2, {_groovy(f2_row["pop1"])}, {_groovy(f2_row["pop2"])})
    POPMIX_F3(f2, {", ".join(_groovy(f3_row[f"pop{i}"]) for i in (1, 2, 3))})
    POPMIX_F4(f2, {", ".join(_groovy(f4_row[f"pop{i}"]) for i in (1, 2, 3, 4))})
    POPMIX_QPWAVE(f2, {_groovy(qpwave["left"])}, {_groovy(qpwave["right"])})
    POPMIX_QPADM(f2, {_groovy(qpadm["target"])}, {_groovy(qpadm["left"])}, {_groovy(qpadm["right"])})
    POPMIX_QPGRAPH(f2.map {{ meta, dir -> [ meta, dir, file('{graph_file}') ] }})

    small = [ [ id:'small' ], {_files(SMALL)} ]
    POPMIX_EXTRACTF2_SMALL(small, {_groovy(SMALL_POPS)})
    POPMIX_FINDGRAPHS(POPMIX_EXTRACTF2_SMALL.out.f2_blocks, {_groovy(SMALL_POPS)})
    POPMIX_DATES(small, 'Simulation', 'French', 'Yoruba')
    POPMIX_CONVERTF([ [ id:'small' ], [ {_files(SMALL)} ] ], 'plink')
    POPMIX_MERGE([ [ id:'a' ], [ {_files(FIXTURES / 'merge_a')} ] ], [ [ id:'b' ], [ {_files(FIXTURES / 'merge_b')} ] ],
                 'eigenstrat')
}}
""")
    (work / "nextflow.config").write_text(f"""
params.outdir = '{work / "results"}'
process {{
    executor = 'local'
    cpus = 4
    publishDir = [ path: {{ "${{params.outdir}}/${{task.process.tokenize(':')[-1].toLowerCase()}}" }}, mode: 'copy' ]
    withName: 'POPMIX_FINDGRAPHS' {{ ext.args = '--nadmix 1 --method exhaustive --top 3' }}
}}
""")
    env = dict(os.environ, PATH=f"{Path(shutil.which('popmix')).parent}:{os.environ['PATH']}", NXF_ANSI_LOG="false")
    result = subprocess.run(["nextflow", "run", "main.nf", "-c", "nextflow.config"], cwd=work, env=env,
                            capture_output=True, text=True, timeout=3600)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]

    def outputs(process: str) -> Path:
        return work / "results" / f"popmix_{process}"

    return outputs


def _json(run, process: str) -> dict:
    (path,) = run(process).glob("*.json")
    return json.loads(path.read_text())


def test_f2(run):
    row, out = _ref("f2.json")["f2"][0], _json(run, "f2")
    assert out["estimate"] == pytest.approx(row["est"], rel=1e-10, abs=1e-15)
    assert out["std_err"] == pytest.approx(row["se"], rel=1e-6)


@pytest.mark.parametrize("stat,n", [("f3", 3), ("f4", 4)])
def test_f3_f4(run, stat, n):
    row, out = _ref("f3_f4.json")[stat][0], _json(run, stat)
    assert out["estimate"] == pytest.approx(row["est"], rel=1e-10, abs=1e-15)
    assert out["std_err"] == pytest.approx(row["se"], rel=1e-6)


def test_qpwave(run):
    m, out = _ref("qpwave.json")["models"][0], _json(run, "qpwave")
    ours = {t["rank"]: t for t in out["rank_tests"]}
    for row in m["rankdrop"]:
        assert ours[row["f4rank"]]["chi_sq"] == pytest.approx(row["chisq"], rel=1e-6)
        assert ours[row["f4rank"]]["p_value"] == pytest.approx(row["p"], rel=1e-4)


def test_qpadm(run):
    m = next(m for m in _ref("qpadm.json")["models"] if not m["constrained"])
    out = _json(run, "qpadm")
    for row in m["weights"]:
        assert out["weights"][row["left"]] == pytest.approx(row["weight"], rel=1e-6)
        assert out["std_errors"][row["left"]] == pytest.approx(row["se"], rel=1e-6)


def test_qpgraph(run):
    graph = next(g for g in _ref("qpgraph.json")["graphs"] if g["name"] == "admix1")
    out = _json(run, "qpgraph")
    assert out["pops"] == graph["pops"]
    assert out["score"] == pytest.approx(graph["score"], rel=1e-10)


def test_findgraphs(run):
    out = _json(run, "findgraphs")
    f2 = pm.extract_f2(str(SMALL), SMALL_POPS)
    graphs, n = pm.find_graphs(f2, SMALL_POPS, nadmix=1, method="exhaustive")
    assert out["n_evaluated"] == n == 84
    assert [g["score"] for g in out["graphs"]] == [g.score for g in graphs[:3]]
    best = (run("findgraphs") / "small_best_graph.txt").read_text().splitlines()
    assert [tuple(line.split("\t")) for line in best if not line.startswith("#")] == graphs[0].edges


def test_dates(run):
    out = _json(run, "dates")
    direct = subprocess.run(["popmix", "dates", "--input", str(SMALL), "--test-pop", "Simulation", "--source1",
                             "French", "--source2", "Yoruba", "--format", "json"],
                            capture_output=True, text=True, check=True)
    assert out == json.loads(direct.stdout)


def test_convertf(run):
    files = sorted(p.name for p in run("convertf").iterdir())
    assert files == ["small_plink.bed", "small_plink.bim", "small_plink.fam"]
    converted = pm.extract_f2(str(run("convertf") / "small_plink"), ["French", "Yoruba"])
    original = pm.extract_f2(str(SMALL), ["French", "Yoruba"])
    assert converted.f2("French", "Yoruba") == original.f2("French", "Yoruba")


def test_merge(run):
    out = run("merge")
    inds = (out / "a_merged.ind").read_text().split()
    assert [w for w in inds if w.startswith("IND")] == ["IND1", "IND2", "IND3"]


SUBWORKFLOW = MODULES.parents[2] / "subworkflows" / "nf-core" / "geno_fstats_admixture_popmix" / "main"
FAMILY = REPO / "tests" / "golden" / "admixtools2" / "v2.0.10" / "family_4pop"


def _f4_reference(pops):
    """admixtools2's f4 for pops, from whichever ordering it lists (f4 changes sign with the pairs)."""
    a, b, c, d = pops
    signs = {(a, b, c, d): 1, (b, a, d, c): 1, (c, d, a, b): 1, (d, c, b, a): 1,
             (b, a, c, d): -1, (a, b, d, c): -1, (d, c, a, b): -1, (c, d, b, a): -1}
    for row in json.loads((FAMILY / "f3_f4.json").read_text())["f4"]:
        key = (row["pop1"], row["pop2"], row["pop3"], row["pop4"])
        if key in signs:
            return signs[key] * row["est"], row["se"]
    raise KeyError(pops)


@pytest.fixture(scope="module")
def subworkflow(tmp_path_factory):
    """Run GENO_FSTATS_ADMIXTURE_POPMIX on the small data; return its published outputs."""
    if not SUBWORKFLOW.with_suffix(".nf").exists():
        pytest.skip(f"no subworkflow at {SUBWORKFLOW}")
    work = tmp_path_factory.mktemp("subworkflow")
    (work / "main.nf").write_text(f"""
include {{ GENO_FSTATS_ADMIXTURE_POPMIX }} from '{SUBWORKFLOW}'

workflow {{
    GENO_FSTATS_ADMIXTURE_POPMIX(
        Channel.of([ [ id:'small' ], {_files(SMALL)} ]),
        'Simulation', ['French', 'Yoruba'], ['Simulation.1', 'Yoruba']
    )
}}
""")
    (work / "nextflow.config").write_text(f"""
params.outdir = '{work / "results"}'
process {{
    executor = 'local'
    cpus = 4
    publishDir = [ path: {{ "${{params.outdir}}/${{task.process.tokenize(':')[-1].toLowerCase()}}" }}, mode: 'copy' ]
}}
""")
    env = dict(os.environ, PATH=f"{Path(shutil.which('popmix')).parent}:{os.environ['PATH']}", NXF_ANSI_LOG="false")
    result = subprocess.run(["nextflow", "run", "main.nf", "-c", "nextflow.config"], cwd=work, env=env,
                            capture_output=True, text=True, timeout=3600)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]
    return lambda process: work / "results" / f"popmix_{process}"


def test_subworkflow_f3_and_f4(subworkflow):
    f3 = _json(subworkflow, "f3")
    ref = next(r for r in json.loads((FAMILY / "f3_f4.json").read_text())["f3"]
               if (r["pop1"], r["pop2"], r["pop3"]) == ("Simulation", "French", "Yoruba"))
    assert f3["estimate"] == pytest.approx(ref["est"], rel=1e-10)
    assert f3["std_err"] == pytest.approx(ref["se"], rel=1e-6)
    est, se = _f4_reference(("Simulation.1", "Simulation", "French", "Yoruba"))
    f4 = _json(subworkflow, "f4")
    assert f4["estimate"] == pytest.approx(est, rel=1e-10)
    assert f4["std_err"] == pytest.approx(se, rel=1e-6)


def test_subworkflow_runs_every_step(subworkflow):
    for process in ["extractf2", "f2", "f3", "f4", "qpadm", "qpwave", "dates"]:
        assert any(subworkflow(process).iterdir()), process
    assert set(_json(subworkflow, "qpadm")["weights"]) == {"French", "Yoruba"}
