"""f-statistic scans: combinations, matrices, clustering order, heatmaps and CLI.

The statistics themselves are checked against admixtools2 in
tests/golden/test_dstat_admixtools2.py; these tests check what is built on
them, on the 54-sample DATES data.
"""

import io
import json
import os
import subprocess

import numpy as np
import pandas as pd
import popmix as pm
import pytest

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def cli(*args: str) -> str:
    return subprocess.run(["popmix", *args], check=True, capture_output=True, text=True).stdout


class TestCombinations:
    def test_all_ordered_triples_without_unique_only(self):
        combs = pm.fstat_combinations(["A", "B"], fnum=3, unique_only=False)
        assert combs[:3] == [("A", "A", "A"), ("A", "A", "B"), ("A", "B", "A")]
        assert len(combs) == 8

    def test_pairs_of_pairs_without_unique_only(self):
        combs = pm.fstat_combinations(["A", "B", "C"], fnum=4, unique_only=False)
        assert combs[0] == ("A", "B", "A", "B")
        assert combs[1] == ("A", "B", "A", "C")
        assert len(combs) == 9

    def test_grid_drops_repeated_combinations(self):
        assert pm.fstat_combinations(["A", "A"], "B", "C", "D", fnum=4) == [("A", "B", "C", "D")]

    def test_too_few_populations_is_an_error(self):
        with pytest.raises(ValueError, match="at least 4"):
            pm.fstat_combinations(["A", "B", "C"], fnum=4)


class TestMatrices:
    def test_outgroup_f3_matrix_is_symmetric_and_matches_qp3pop(self):
        others = ["French", "Simulation", "Simulation.1"]
        m = pm.outgroup_f3_matrix(SMALL, "Yoruba", others)
        assert list(m.index) == others and list(m.columns) == others
        np.testing.assert_array_equal(m.values, m.values.T)
        rows = pm.qp3pop(SMALL, "Yoruba", ["French"], ["Simulation", "French"], outgroupmode=True)
        assert m.loc["French", "Simulation"] == rows["est"][0]
        assert m.loc["French", "French"] == rows["est"][1]

    def test_outgroup_f3_z_matrix(self):
        m = pm.outgroup_f3_matrix(SMALL, "Yoruba", ["French", "Simulation"], value="z")
        rows = pm.qp3pop(SMALL, "Yoruba", "French", "Simulation", outgroupmode=True)
        assert m.loc["Simulation", "French"] == rows["z"][0]

    def test_dstat_matrix_is_antisymmetric_and_matches_qpdstat(self):
        others = ["French", "Simulation", "Simulation.1"]
        m = pm.dstat_matrix(SMALL, others, pop3="Yoruba", pop4="French")
        np.testing.assert_array_equal(m.values, -m.values.T)
        assert np.isnan(np.diag(m.values)).all()
        rows = pm.qpdstat(SMALL, "Simulation", "Simulation.1", "Yoruba", "French", f4mode=False)
        assert m.loc["Simulation", "Simulation.1"] == rows["est"][0]

    def test_dstat_matrix_over_the_second_pair(self):
        m = pm.dstat_matrix(SMALL, ["Simulation", "Simulation.1"], pop1="Yoruba", pop2="French", value="z")
        rows = pm.qpdstat(SMALL, "Yoruba", "French", "Simulation", "Simulation.1", f4mode=False)
        assert m.loc["Simulation", "Simulation.1"] == rows["z"][0]

    def test_dstat_matrix_needs_exactly_one_fixed_pair(self):
        with pytest.raises(ValueError, match="pop1 and pop2, or pop3 and pop4"):
            pm.dstat_matrix(SMALL, POPS, pop1="Yoruba")

    def test_f4_matrix_from_f2_blocks(self):
        f2 = pm.extract_f2(SMALL, POPS)
        m = pm.dstat_matrix(f2, ["Simulation", "Simulation.1"], pop3="Yoruba", pop4="French", f4mode=True)
        assert m.loc["Simulation", "Simulation.1"] == f2.f4("Simulation", "Simulation.1", "Yoruba", "French")["estimate"]


class TestClusterOrder:
    def test_similar_rows_end_up_adjacent(self):
        names = ["a1", "b1", "a2", "b2", "a3"]
        a, b = [1.0, 0.0, 1.0, 0.0, 1.0], [0.0, 1.0, 0.0, 1.0, 0.0]
        rows = [a, b, [x + 0.01 for x in a], [x + 0.01 for x in b], [x - 0.01 for x in a]]
        order = pm.cluster_order(pd.DataFrame(rows, index=names, columns=names))
        groups = "".join(n[0] for n in order)
        assert groups in ("aaabb", "bbaaa")

    def test_order_is_deterministic_and_complete(self):
        m = pm.outgroup_f3_matrix(SMALL, "Yoruba", ["French", "Simulation", "Simulation.1"])
        assert sorted(pm.cluster_order(m)) == sorted(m.index)
        assert pm.cluster_order(m) == pm.cluster_order(m)


class TestHeatmap:
    def test_heatmap_writes_an_image(self, tmp_path):
        pytest.importorskip("matplotlib")
        from popmix import plot

        m = pm.outgroup_f3_matrix(SMALL, "Yoruba", ["French", "Simulation", "Simulation.1"])
        out = tmp_path / "f3.png"
        plot.heatmap(m, title="f3(Yoruba; X, Y)", path=out)
        assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


class TestCli:
    def test_qpdstat_tsv_matches_the_api(self):
        out = cli("qpdstat", "--input", SMALL, "--pop1", ",".join(POPS), "--d-mode", "--format", "tsv")
        ours = pd.read_csv(io.StringIO(out), sep="\t")
        api = pm.qpdstat(SMALL, POPS, f4mode=False)
        assert ours[["pop1", "pop2", "pop3", "pop4"]].values.tolist() == api[["pop1", "pop2", "pop3", "pop4"]].values.tolist()
        # TSV has 12 significant digits
        np.testing.assert_allclose(ours["est"], api["est"], rtol=1e-11)
        assert ours["n"].tolist() == api["n"].tolist()

    def test_qpdstat_from_saved_f2_blocks(self, tmp_path):
        cli("extract-f2", "--input", SMALL, "--populations", ",".join(POPS), "--output", str(tmp_path / "f2"))
        out = json.loads(cli("qpdstat", "--f2-dir", str(tmp_path / "f2"), "--pop1", "Simulation",
                             "--pop2", "Simulation.1", "--pop3", "Yoruba", "--pop4", "French", "--format", "json"))
        f2 = pm.read_f2(str(tmp_path / "f2"))
        assert out[0]["est"] == pytest.approx(f2.f4("Simulation", "Simulation.1", "Yoruba", "French")["estimate"],
                                              rel=1e-12)

    def test_d_mode_needs_genotypes(self, tmp_path):
        cli("extract-f2", "--input", SMALL, "--populations", ",".join(POPS), "--output", str(tmp_path / "f2"))
        r = subprocess.run(["popmix", "qpdstat", "--f2-dir", str(tmp_path / "f2"), "--pop1", ",".join(POPS),
                            "--d-mode"], capture_output=True, text=True)
        assert r.returncode != 0
        assert "need genotypes" in r.stderr

    def test_qp3pop_outgroup_mode(self):
        out = cli("qp3pop", "--input", SMALL, "--pop1", "Yoruba", "--pop2", "French", "--pop3", "Simulation",
                  "--outgroup-mode", "--format", "json")
        api = pm.qp3pop(SMALL, "Yoruba", "French", "Simulation", outgroupmode=True)
        assert json.loads(out)[0]["est"] == pytest.approx(api["est"][0], rel=1e-15)

    def test_outgroup_f3_matrix_and_heatmap(self, tmp_path):
        pytest.importorskip("matplotlib")
        png = tmp_path / "f3.png"
        out = cli("outgroup-f3", "--input", SMALL, "--outgroup", "Yoruba", "--populations",
                  "French,Simulation,Simulation.1", "--heatmap", str(png), "--order", "cluster")
        m = pd.read_csv(io.StringIO(out), sep="\t", index_col=0)
        api = pm.outgroup_f3_matrix(SMALL, "Yoruba", ["French", "Simulation", "Simulation.1"])
        order = pm.cluster_order(api)
        assert list(m.index) == order
        np.testing.assert_allclose(m.values, api.loc[order, order].values, rtol=1e-11)
        assert png.read_bytes()[:4] == b"\x89PNG"

    def test_dstat_matrix_command(self):
        out = cli("dstat-matrix", "--input", SMALL, "--populations", "French,Simulation,Simulation.1",
                  "--pop3", "Yoruba", "--pop4", "French", "--value", "z")
        m = pd.read_csv(io.StringIO(out), sep="\t", index_col=0)
        api = pm.dstat_matrix(SMALL, ["French", "Simulation", "Simulation.1"], pop3="Yoruba", pop4="French",
                              value="z")
        np.testing.assert_allclose(m.values, api.values, rtol=1e-11)
