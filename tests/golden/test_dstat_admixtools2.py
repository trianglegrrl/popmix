"""Golden master tests: f4-, D- and f3-statistics from genotype files vs admixtools2.

The reference (tests/golden/admixtools2/v<version>/<config>/dstat.json) is
admixtools2's qpdstat() and qp3pop() given a genotype prefix, written by
generate_dstat.R. These use every autosomal SNP, blocks formed over all of
them, and each statistic's own SNPs (allsnps = TRUE), so they differ from the
same statistics computed from f2 blocks.

Stages, so a mismatch points at its cause:

- blocks: per-block numerator, denominator and SNP count for a few
  combinations (f4blockdat_from_geno, f3blockdat_from_geno), rtol 1e-10 and
  counts exact
- combinations: the population combinations qpdstat and qp3pop generate, in
  their order
- qpdstat_f4, qpdstat_d, qpdstat_grid_d, qp3pop, qp3pop_outgroup: estimates
  (rtol 1e-10), standard errors and z (rtol 1e-6), p-values (rtol 1e-4) and
  SNP counts (exact)

mixed_ploidy relabels AADR samples so one population mixes pseudo-haploid and
diploid samples, where frequencies (every sample counted as diploid) and the
f3 correction (pseudo-haploid samples counted once) differ.
"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import popmix as pm

from .reference_paths import local, missing

REF_ROOT = Path(__file__).parent / "admixtools2" / "v2.0.10"
# The AADR configs read the 5 GB genotype file several times: --run-slow
CONFIGS = ["family_4pop"] + [pytest.param(c, marks=pytest.mark.slow)
                             for c in ["dg_6pop", "ancient_5pop", "mixed_ploidy"]]

F_RTOL = 1e-10
F_ATOL = 1e-15
SE_RTOL = 1e-6
P_RTOL = 1e-4


def _reference(name: str) -> dict:
    path = REF_ROOT / name / "dstat.json"
    if not path.exists():
        pytest.skip(f"no reference: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


_prefixes: dict = {}


def _prefix(name: str, ref: dict, tmp_root: Path) -> str:
    """The genotype prefix, relabeled as the reference was (symlinked .geno/.snp)."""
    if missing(f"{local(ref['prefix'])}.geno"):
        pytest.skip(f"genotype data not available: {local(ref['prefix'])}")
    if not ref.get("relabel"):
        return local(ref["prefix"])
    if name not in _prefixes:
        out = tmp_root / name / Path(local(ref["prefix"])).name
        out.parent.mkdir(parents=True)
        for ext in (".geno", ".snp"):
            os.symlink(os.path.realpath(local(ref["prefix"]) + ext), str(out) + ext)
        relabel = {old: new for new, olds in ref["relabel"].items() for old in olds}
        lines = []
        for line in Path(local(ref["prefix"]) + ".ind").read_text(encoding="utf-8").splitlines():
            ind, sex, pop = line.split()
            lines.append(f"{ind}\t{sex}\t{relabel.get(pop, pop)}")
        Path(str(out) + ".ind").write_text("\n".join(lines) + "\n", encoding="utf-8")
        _prefixes[name] = str(out)
    return _prefixes[name]


@pytest.fixture(scope="module")
def tmp_root(tmp_path_factory) -> Path:
    return tmp_path_factory.mktemp("relabeled")


def _assert_stats(ours: pd.DataFrame, rows: list, pop_cols: list):
    theirs = pd.DataFrame(rows)
    assert ours[pop_cols].values.tolist() == theirs[pop_cols].values.tolist()
    # NaN (JSON null) where z is 0/0, as for f(A, A; C, D); p is 0 where it underflows
    num = {k: theirs[k].astype(float) for k in ["est", "se", "z", "p"]}
    np.testing.assert_allclose(ours["est"], num["est"], rtol=F_RTOL, atol=F_ATOL)
    np.testing.assert_allclose(ours["se"], num["se"], rtol=SE_RTOL, atol=0)
    np.testing.assert_allclose(ours["z"], num["z"], rtol=SE_RTOL, atol=0)
    np.testing.assert_allclose(ours["p"], num["p"], rtol=P_RTOL, atol=1e-300)
    # admixtools2 gives NA for n where nothing can be estimated; popmix gives 0
    assert ours["n"].tolist() == theirs["n"].fillna(0).astype(int).tolist()
    assert np.isnan(num["est"][ours["n"] == 0]).all()


@pytest.mark.parametrize("name", CONFIGS)
def test_ploidy(name, tmp_root):
    ref = _reference(name)
    prefix = _prefix(name, ref, tmp_root)
    assert pm.sample_ploidy(prefix, ref["pops"]) == ref["ploidy"]


@pytest.mark.parametrize("name", CONFIGS)
def test_f4_blocks(name, tmp_root):
    ref = _reference(name)
    theirs = pd.DataFrame(ref["f4_blocks"])
    combs = theirs[["pop1", "pop2", "pop3", "pop4"]].drop_duplicates().values.tolist()
    ours = pm.f4blockdat_from_geno(_prefix(name, ref, tmp_root), combs, f4mode=False)
    cols = ["pop1", "pop2", "pop3", "pop4", "block"]
    assert ours[cols].values.tolist() == theirs[cols].values.tolist()
    assert ours["n"].tolist() == theirs["n"].tolist()
    assert ours["length"].tolist() == theirs["length"].tolist()
    np.testing.assert_allclose(ours["est"], theirs["est"].astype(float), rtol=F_RTOL, atol=F_ATOL)
    np.testing.assert_allclose(ours["den"], theirs["den"].astype(float), rtol=F_RTOL, atol=F_ATOL)


@pytest.mark.parametrize("name", CONFIGS)
def test_f3_blocks(name, tmp_root):
    ref = _reference(name)
    theirs = pd.DataFrame(ref["f3_blocks"])
    combs = theirs[["pop1", "pop2", "pop3"]].drop_duplicates().values.tolist()
    ours = pm.f3blockdat_from_geno(_prefix(name, ref, tmp_root), combs)
    cols = ["pop1", "pop2", "pop3", "block"]
    assert ours[cols].values.tolist() == theirs[cols].values.tolist()
    assert ours["n"].tolist() == theirs["n"].tolist()
    for col in ["numer", "denom", "est"]:
        np.testing.assert_allclose(ours[col], theirs[col].astype(float), rtol=F_RTOL, atol=F_ATOL)


@pytest.mark.parametrize("name", CONFIGS)
def test_combinations(name):
    ref = _reference(name)
    pops, g = ref["pops"], ref["grid"]
    four = ["pop1", "pop2", "pop3", "pop4"]
    assert pm.fstat_combinations(pops, fnum=4) == [tuple(r) for r in pd.DataFrame(ref["qpdstat_d"])[four].values]
    assert pm.fstat_combinations(pops, fnum=3) == [
        tuple(r) for r in pd.DataFrame(ref["qp3pop"])[four[:3]].values
    ]
    grid = pm.fstat_combinations(g["pop1"], g["pop2"], g["pop3"], g["pop4"], fnum=4)
    assert grid == [tuple(r) for r in pd.DataFrame(ref["qpdstat_grid_d"])[four].values]


@pytest.mark.parametrize("name", CONFIGS)
def test_qpdstat_f4(name, tmp_root):
    ref = _reference(name)
    ours = pm.qpdstat(_prefix(name, ref, tmp_root), ref["pops"])
    _assert_stats(ours, ref["qpdstat_f4"], ["pop1", "pop2", "pop3", "pop4"])


@pytest.mark.parametrize("name", CONFIGS)
def test_qpdstat_d(name, tmp_root):
    ref = _reference(name)
    ours = pm.qpdstat(_prefix(name, ref, tmp_root), ref["pops"], f4mode=False)
    _assert_stats(ours, ref["qpdstat_d"], ["pop1", "pop2", "pop3", "pop4"])


@pytest.mark.parametrize("name", CONFIGS)
def test_qpdstat_grid_d(name, tmp_root):
    ref = _reference(name)
    g = ref["grid"]
    ours = pm.qpdstat(_prefix(name, ref, tmp_root), g["pop1"], g["pop2"], g["pop3"], g["pop4"], f4mode=False)
    _assert_stats(ours, ref["qpdstat_grid_d"], ["pop1", "pop2", "pop3", "pop4"])


@pytest.mark.parametrize("name", CONFIGS)
def test_qp3pop(name, tmp_root):
    ref = _reference(name)
    ours = pm.qp3pop(_prefix(name, ref, tmp_root), ref["pops"])
    _assert_stats(ours, ref["qp3pop"], ["pop1", "pop2", "pop3"])


@pytest.mark.parametrize("name", CONFIGS)
def test_qp3pop_outgroup(name, tmp_root):
    ref = _reference(name)
    ours = pm.qp3pop(_prefix(name, ref, tmp_root), ref["pops"], outgroupmode=True)
    _assert_stats(ours, ref["qp3pop_outgroup"], ["pop1", "pop2", "pop3"])


@pytest.mark.parametrize("name", ["family_4pop", pytest.param("dg_6pop", marks=pytest.mark.slow)])
def test_scans_from_f2_blocks(name):
    """qpdstat and qp3pop on f2 blocks: admixtools2's f4() and f3() on its f2 blocks (f3_f4.json)."""
    ref = json.loads((REF_ROOT / name / "f2.json").read_text(encoding="utf-8"))
    if missing(f"{local(ref['prefix'])}.geno"):
        pytest.skip(f"genotype data not available: {local(ref['prefix'])}")
    f2 = pm.extract_f2(local(ref["prefix"]), ref["pops"])
    stats = json.loads((REF_ROOT / name / "f3_f4.json").read_text(encoding="utf-8"))
    # admixtools2 returns qp3pop() on f2 blocks sorted by population (dplyr
    # summarize); popmix keeps the combination order, as both do from genotypes.
    f3_cols = ["pop1", "pop2", "pop3"]
    for ours, rows, cols in [(pm.qpdstat(f2), stats["f4"], ["pop1", "pop2", "pop3", "pop4"]),
                             (pm.qp3pop(f2).sort_values(f3_cols, ignore_index=True), stats["f3"], f3_cols)]:
        theirs = pd.DataFrame(rows)
        assert ours[cols].values.tolist() == theirs[cols].values.tolist()
        np.testing.assert_allclose(ours["est"], theirs["est"], rtol=F_RTOL, atol=F_ATOL)
        np.testing.assert_allclose(ours["se"], theirs["se"], rtol=SE_RTOL, atol=0)
