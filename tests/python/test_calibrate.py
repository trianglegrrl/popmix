"""Calibration by simulation (TASKS.md H7).

The replicates are drawn here with numpy, a Balding-Nichols tree
((A, B), (C, D)) with outgroups E and F, with or without gene flow from A
into C, written as
EIGENSTRAT; users make theirs with msprime, which popmix never imports.
On the tree f4(A, B; C, D) is 0, and the rejection rate is the false-
positive rate; with gene flow it is the power.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import popmix as pm
from popmix.calibrate import wilson_interval

N_SNPS, N_CHROM, PER_POP = 3000, 10, 10


def drift(p: np.ndarray, fst: float, rng: np.random.Generator) -> np.ndarray:
    a = p * (1 - fst) / fst
    return rng.beta(a, (1 - p) * (1 - fst) / fst)


def replicate(path: Path, seed: int, admixture: float) -> str:
    rng = np.random.default_rng(seed)
    p0 = rng.uniform(0.1, 0.9, N_SNPS)
    ab, cd = drift(p0, 0.05, rng), drift(p0, 0.05, rng)
    a, b, c, d = drift(ab, 0.05, rng), drift(ab, 0.05, rng), drift(cd, 0.05, rng), drift(cd, 0.05, rng)
    e, f = drift(p0, 0.1, rng), drift(p0, 0.1, rng)
    c = (1 - admixture) * c + admixture * a
    geno = np.concatenate([rng.binomial(2, x[:, None], (N_SNPS, PER_POP)) for x in (a, b, c, d, e, f)], axis=1)
    prefix = path / f"rep{seed}"
    per = N_SNPS // N_CHROM
    Path(f"{prefix}.snp").write_text("".join(
        f"s{k} {k // per + 1} {(k % per) * 0.5 / per:.6f} {(k % per) * 1000 + 1} A G\n" for k in range(N_SNPS)))
    Path(f"{prefix}.geno").write_text("".join("".join(map(str, row)) + "\n" for row in geno))
    Path(f"{prefix}.ind").write_text("".join(f"{p}{i} U {p}\n" for p in "ABCDEF" for i in range(PER_POP)))
    return str(prefix)


@pytest.fixture(scope="module")
def null(tmp_path_factory) -> list[str]:
    d = tmp_path_factory.mktemp("null")
    return [replicate(d, s, 0.0) for s in range(20)]


@pytest.fixture(scope="module")
def flow(tmp_path_factory) -> list[str]:
    d = tmp_path_factory.mktemp("flow")
    return [replicate(d, 100 + s, 0.3) for s in range(10)]


def test_the_wilson_interval():
    lo, hi = wilson_interval(5, 100)
    assert (round(lo, 4), round(hi, 4)) == (0.0215, 0.1118)
    assert wilson_interval(0, 0) == (0.0, 1.0)


def test_a_tree_rejects_at_about_alpha(null):
    df = pm.calibrate_f4(null, "A", "B", "C", "D", block_size=0.02)
    assert df.attrs["n"] == 20
    assert df.attrs["ci_low"] <= 0.05 <= df.attrs["ci_high"]
    assert df["rejected"].tolist() == (df["p"] <= 0.05).tolist()


def test_gene_flow_is_found(flow):
    df = pm.calibrate_f4(flow, "A", "B", "C", "D", block_size=0.02)
    assert df.attrs["rate"] >= 0.8


def test_qpadm_on_the_true_sources(flow):
    # C is a mix of the C-D branch and A; D stands in for its unadmixed side
    df = pm.calibrate_qpadm(flow, ["A", "D"], ["E", "F", "B"], "C", block_size=0.02)
    assert len(df) == 10
    assert set(df.columns) == {"replicate", "chisq", "dof", "p", "feasible", "weights", "rejected"}


def test_the_cli_summarizes(flow, tmp_path):
    r = subprocess.run(["popmix", "calibrate-f4", *flow, "--pop1", "A", "--pop2", "B", "--pop3", "C", "--pop4", "D",
                        "--block-size", "0.02", "-o", str(tmp_path / "reps.tsv")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    a = pm.calibrate_f4(flow, "A", "B", "C", "D", block_size=0.02).attrs
    assert (f"10 replicates, {a['rejected']} rejected at alpha 0.05: rate {a['rate']:.3f} "
            f"(95% CI {a['ci_low']:.3f}-{a['ci_high']:.3f})") in r.stdout
    assert len((tmp_path / "reps.tsv").read_text().splitlines()) == 11


def test_popmix_never_imports_msprime():
    code = "import popmix, popmix.calibrate, popmix.cli, sys; print('msprime' in sys.modules)"
    assert subprocess.run([sys.executable, "-c", code], capture_output=True, text=True).stdout.strip() == "False"
