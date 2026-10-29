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
        f"s{k} {k // per + 1} {(k % per) * 0.5 / per:.6f} {(k % per) * 1000 + 1} A G\n" for k in range(N_SNPS)), encoding="utf-8")
    Path(f"{prefix}.geno").write_text("".join("".join(map(str, row)) + "\n" for row in geno), encoding="utf-8")
    Path(f"{prefix}.ind").write_text("".join(f"{p}{i} U {p}\n" for p in "ABCDEF" for i in range(PER_POP)), encoding="utf-8")
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
    assert len((tmp_path / "reps.tsv").read_text(encoding="utf-8").splitlines()) == 11


def test_popmix_never_imports_msprime():
    code = "import popmix, popmix.calibrate, popmix.cli, sys; print('msprime' in sys.modules)"
    assert subprocess.run([sys.executable, "-c", code], capture_output=True, text=True).stdout.strip() == "False"


# What the tests above left unchecked (TASKS.md K6: each kills a mutant that
# survived them).

@pytest.mark.parametrize("n", [9, 26, 30])
def test_the_wilson_interval_of_all_rejected_ends_at_1(n):
    assert wilson_interval(n, n)[1] == 1.0


@pytest.mark.parametrize("n", [10, 21, 22])
def test_the_wilson_interval_of_none_rejected_starts_at_0(n):
    assert wilson_interval(0, n)[0] == 0.0


def test_the_rate_is_the_share_rejected(null):
    a = pm.calibrate_f4(null, "A", "B", "C", "D", block_size=0.02, alpha=0.5).attrs
    assert a["rejected"] > 1
    assert a["rate"] == a["rejected"] / 20


@pytest.mark.parametrize("call", [lambda: pm.calibrate_f4([], "A", "B", "C", "D"),
                                  lambda: pm.calibrate_qpadm([], ["A", "D"], ["E", "F", "B"], "C")])
def test_no_replicates_is_an_error(call):
    with pytest.raises(ValueError, match="^no replicates to calibrate on$"):
        call()


def test_the_f4_columns(null):
    df = pm.calibrate_f4(null[:2], "A", "B", "C", "D", block_size=0.02)
    assert df.columns.tolist() == ["replicate", "est", "se", "z", "p", "rejected"]


def test_an_f4_p_at_alpha_rejects(null):
    p = pm.calibrate_f4(null[:1], "A", "B", "C", "D", block_size=0.02)["p"][0]
    assert pm.calibrate_f4(null[:1], "A", "B", "C", "D", block_size=0.02, alpha=p)["rejected"].tolist() == [True]


@pytest.fixture(scope="module")
def flow_f2(flow) -> list:
    return [pm.extract_f2(rep, list("ABCDEF"), block_size=0.02) for rep in flow[:2]]


def test_qpadm_calibration_defaults_to_alpha_0_05(flow_f2):
    assert pm.calibrate_qpadm(flow_f2, ["A", "D"], ["E", "F", "B"], "C").attrs["alpha"] == 0.05


def test_qpadm_calibration_is_unconstrained_by_default(flow_f2):
    df = pm.calibrate_qpadm(flow_f2[:1], ["A", "D"], ["E", "F", "B"], "C")
    r = pm.qpadm(flow_f2[0], ["A", "D"], ["E", "F", "B"], "C")
    assert df["weights"][0] == f"{r.weights[0]:.6g},{r.weights[1]:.6g}"


def test_qpadm_calibration_differs_constrained_when_a_weight_is_out_of_range(flow_f2):
    # B from C and D, a weight outside [0, 1] (-0.549) on the first replicate
    r = pm.qpadm(flow_f2[0], ["C", "D"], ["A", "E", "F"], "B")
    assert not r.feasible
    df = pm.calibrate_qpadm(flow_f2[:1], ["C", "D"], ["A", "E", "F"], "B")
    assert df["weights"][0] == f"{r.weights[0]:.6g},{r.weights[1]:.6g}"


def test_a_qpadm_p_at_alpha_rejects(flow_f2):
    p = pm.qpadm(flow_f2[0], ["A", "D"], ["E", "F", "B"], "C").p
    df = pm.calibrate_qpadm(flow_f2[:1], ["A", "D"], ["E", "F", "B"], "C", alpha=p)
    assert df["rejected"].tolist() == [True]
