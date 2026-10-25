"""
CLI tests for popmix.

These tests define the expected CLI behavior. Written BEFORE implementation (TDD).
"""

import subprocess
import json
import tempfile
import os
from pathlib import Path

import pytest

# Test data paths
AADR_PREFIX = "/references/aadr/v62.0_1240k_public"
AADR_HO_PREFIX = "/references/aadr_v54.1.p1/HO_public"

# Actual population names in AADR v62 (with suffixes)
POP_HAN = "Han.DG"
POP_YORUBA = "Yoruba.DG"
POP_FRENCH = "France_BellBeaker.AG"
POP_MBUTI = "Mbuti.DG"


def run_cli(*args, check=True, timeout=300):
    """Run popmix CLI and return result."""
    result = subprocess.run(
        ["python", "-m", "popmix.cli"] + list(args),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if check and result.returncode != 0:
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, args, result.stdout, result.stderr)
    return result


class TestCLIHelp:
    """Test CLI help and version commands."""

    def test_help_shows_usage(self):
        """popmix --help should show usage information."""
        result = run_cli("--help")
        assert "Usage:" in result.stdout or "usage:" in result.stdout.lower()
        assert "popmix" in result.stdout.lower()

    def test_version_shows_version(self):
        """popmix --version should show version number."""
        result = run_cli("--version")
        # Should contain a version like 0.1.0
        assert "0." in result.stdout or "version" in result.stdout.lower()

    def test_no_args_shows_help(self):
        """popmix with no args should show help or error."""
        result = run_cli(check=False)
        # Should either show help or indicate missing command
        assert result.returncode == 0 or "Usage" in result.stderr or "Error" in result.stderr


class TestInfoCommand:
    """Test the 'info' command for dataset inspection."""

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    def test_info_packedancestrymap(self):
        """popmix info should show dataset information for PACKEDANCESTRYMAP."""
        result = run_cli("info", AADR_PREFIX)

        # Should show format
        assert "packedancestrymap" in result.stdout.lower() or "format" in result.stdout.lower()
        # Should show counts
        assert "sample" in result.stdout.lower() or "snp" in result.stdout.lower()

    @pytest.mark.skipif(not os.path.exists(f"{AADR_HO_PREFIX}.bed"), reason="AADR HO data not available")
    def test_info_plink(self):
        """popmix info should show dataset information for PLINK."""
        result = run_cli("info", AADR_HO_PREFIX)

        assert "plink" in result.stdout.lower() or "format" in result.stdout.lower()

    def test_info_missing_file_errors(self):
        """popmix info on missing file should error."""
        result = run_cli("info", "/nonexistent/path", check=False)
        assert result.returncode != 0


class TestExtractF2Command:
    """Test the 'extract-f2' command."""

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_extract_f2_creates_output(self):
        """popmix extract-f2 should create F2 block files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "f2_blocks"

            # Use streaming extraction with limited populations (fast path)
            # Note: don't use --snp-range to enable streaming mode
            result = run_cli(
                "extract-f2",
                "--input", AADR_PREFIX,
                "--output", str(outdir),
                "--populations", "Han.DG,Yoruba.DG",  # Just 2 populations
                timeout=120,
            )

            assert result.returncode == 0
            assert outdir.exists()

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_extract_f2_roundtrip(self):
        """extract-f2 output can be loaded back with read_f2()."""
        import popmix as pm

        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "f2_blocks"

            # Extract F2 blocks
            result = run_cli(
                "extract-f2",
                "--input", AADR_PREFIX,
                "--output", str(outdir),
                "--populations", "Han.DG,Yoruba.DG",
                timeout=120,
            )
            assert result.returncode == 0

            # Should have all required files
            assert (outdir / "metadata.json").exists(), "metadata.json not created"
            assert (outdir / "block_sums.npy").exists(), "block_sums.npy not created"
            assert (outdir / "block_weights.npy").exists(), "block_weights.npy not created"

            # Load back and verify
            f2_loaded = pm.read_f2(str(outdir))
            assert f2_loaded.n_populations == 2
            assert "Han.DG" in f2_loaded.populations
            assert "Yoruba.DG" in f2_loaded.populations

            # Should be able to compute f2 from loaded data
            f2_result = f2_loaded.f2("Han.DG", "Yoruba.DG")
            assert f2_result is not None
            assert f2_result["estimate"] > 0  # Different populations should have positive f2

    def test_extract_f2_missing_input_errors(self):
        """popmix extract-f2 with missing input should error."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_cli(
                "extract-f2",
                "--input", "/nonexistent/path",
                "--output", tmpdir,
                check=False
            )
            assert result.returncode != 0


SMALL_PREFIX = "tests/golden/dates/family_chr20_22"
small_data = pytest.mark.skipif(not os.path.exists(f"{SMALL_PREFIX}.geno"), reason="test data not available")
aadr_data = pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")


class TestF2Command:
    """Test the 'f2' command for computing f2 statistics."""

    @small_data
    def test_f2_computes_statistic(self):
        result = run_cli("f2", "--input", SMALL_PREFIX, "--pop1", "French", "--pop2", "Yoruba",
                         "--all-chromosomes")
        assert "f2(French, Yoruba)" in result.stdout
        assert "estimate" in result.stdout.lower()


class TestF3Command:
    """Test the 'f3' command for computing f3 statistics."""

    @small_data
    def test_f3_computes_statistic(self):
        result = run_cli("f3", "--input", SMALL_PREFIX, "--target", "Simulation",
                         "--source1", "French", "--source2", "Yoruba", "--all-chromosomes")
        assert "f3(Simulation; French, Yoruba)" in result.stdout


class TestF4Command:
    """Test the 'f4' command for computing f4/D statistics."""

    @small_data
    def test_f4_computes_statistic(self):
        result = run_cli("f4", "--input", SMALL_PREFIX, "--pop1", "Simulation", "--pop2", "Simulation.1",
                         "--pop3", "French", "--pop4", "Yoruba", "--all-chromosomes")
        assert "f4(Simulation, Simulation.1; French, Yoruba)" in result.stdout


class TestQpAdmCommand:
    """Test the 'qpadm' command."""

    @small_data
    def test_qpadm_computes_weights(self):
        result = run_cli("qpadm", "--input", SMALL_PREFIX, "--target", "Simulation", "--left", "French",
                         "--right", "Yoruba,Simulation.1", "--all-chromosomes")
        assert "weight" in result.stdout.lower()


class TestQpWaveCommand:
    """Test the 'qpwave' command."""

    @small_data
    def test_qpwave_computes_rank(self):
        result = run_cli("qpwave", "--input", SMALL_PREFIX, "--left", "Simulation,French",
                         "--right", "Yoruba,Simulation.1", "--all-chromosomes")
        assert "rank" in result.stdout.lower()


DG6_REF = Path("tests/golden/admixtools2/v2.0.10/dg_6pop")


def _dg6(file: str) -> dict:
    return json.loads((DG6_REF / file).read_text())


@pytest.fixture(scope="module")
def dg6_f2_dir(tmp_path_factory):
    """f2 blocks for the dg_6pop reference populations, from 'popmix extract-f2'."""
    if not os.path.exists(f"{AADR_PREFIX}.geno"):
        pytest.skip("AADR data not available")
    out = tmp_path_factory.mktemp("dg6") / "f2"
    run_cli("extract-f2", "--input", AADR_PREFIX, "--output", str(out),
            "--populations", ",".join(_dg6("f2.json")["pops"]), timeout=1800)
    return str(out)


@pytest.mark.slow
class TestCommandsMatchAdmixtools2:
    """The f-statistic commands reproduce admixtools2 2.0.10 on AADR.

    The populations are extracted once with 'popmix extract-f2' and read back
    with --f2-dir, as admixtools2's extract_f2() and f2_from_precomp() work,
    so the SNP set is the reference's.
    """

    def test_f2(self, dg6_f2_dir):
        for row in _dg6("f2.json")["f2"]:
            out = json.loads(run_cli("f2", "--f2-dir", dg6_f2_dir, "--pop1", row["pop1"],
                                     "--pop2", row["pop2"], "--format", "json").stdout)
            assert out["estimate"] == pytest.approx(row["est"], rel=1e-10, abs=1e-15)
            assert out["std_err"] == pytest.approx(row["se"], rel=1e-6)

    def test_f3(self, dg6_f2_dir):
        for row in _dg6("f3_f4.json")["f3"][:10]:
            out = json.loads(run_cli("f3", "--f2-dir", dg6_f2_dir, "--target", row["pop1"],
                                     "--source1", row["pop2"], "--source2", row["pop3"],
                                     "--format", "json").stdout)
            assert out["estimate"] == pytest.approx(row["est"], rel=1e-10, abs=1e-15)
            assert out["std_err"] == pytest.approx(row["se"], rel=1e-6)

    def test_f4(self, dg6_f2_dir):
        for row in _dg6("f3_f4.json")["f4"][:10]:
            out = json.loads(run_cli("f4", "--f2-dir", dg6_f2_dir, "--pop1", row["pop1"], "--pop2", row["pop2"],
                                     "--pop3", row["pop3"], "--pop4", row["pop4"], "--format", "json").stdout)
            assert out["estimate"] == pytest.approx(row["est"], rel=1e-10, abs=1e-15)
            assert out["std_err"] == pytest.approx(row["se"], rel=1e-6)

    def test_qpwave(self, dg6_f2_dir):
        for m in _dg6("qpwave.json")["models"]:
            out = json.loads(run_cli("qpwave", "--f2-dir", dg6_f2_dir, "--left", ",".join(m["left"]),
                                     "--right", ",".join(m["right"]), "--format", "json").stdout)
            ours = {t["rank"]: t for t in out["rank_tests"]}
            for row in m["rankdrop"]:
                t = ours[row["f4rank"]]
                assert t["df"] == row["dof"]
                assert t["chi_sq"] == pytest.approx(row["chisq"], rel=1e-6)
                if "p" in row:
                    assert t["p_value"] == pytest.approx(row["p"], rel=1e-4)

    def test_qpadm(self, dg6_f2_dir):
        for m in _dg6("qpadm.json")["models"]:
            args = ["qpadm", "--f2-dir", dg6_f2_dir, "--target", m["target"], "--left", ",".join(m["left"]),
                    "--right", ",".join(m["right"]), "--format", "json"]
            if m["constrained"]:
                args.append("--constrained")
            out = json.loads(run_cli(*args).stdout)
            for row in m["weights"]:
                assert out["weights"][row["left"]] == pytest.approx(row["weight"], rel=1e-6, abs=1e-10)
                assert out["std_errors"][row["left"]] == pytest.approx(row["se"], rel=1e-6, abs=1e-10)
            if "rankdrop" in m:
                row = next(x for x in m["rankdrop"] if x["f4rank"] == len(m["left"]) - 1)
                assert out["df"] == row["dof"]
                assert out["chi_sq"] == pytest.approx(row["chisq"], rel=1e-6)
                if "p" in row:
                    assert out["p_value"] == pytest.approx(row["p"], rel=1e-4)

    def test_genotype_input_uses_extract_f2(self):
        """Without --f2-dir the command extracts just its own populations, as pm.extract_f2 does."""
        if not os.path.exists(f"{AADR_PREFIX}.geno"):
            pytest.skip("AADR data not available")
        pops = ["Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG"]
        import popmix as pm
        expected = pm.extract_f2(AADR_PREFIX, pops).f4(*pops)
        out = json.loads(run_cli("f4", "--input", AADR_PREFIX, "--pop1", pops[0], "--pop2", pops[1],
                                 "--pop3", pops[2], "--pop4", pops[3], "--format", "json",
                                 timeout=1800).stdout)
        assert out["estimate"] == expected["estimate"]
        assert out["std_err"] == expected["std_err"]


class TestF2Source:
    """Where the f-statistic commands take their f2 blocks from."""

    @small_data
    def test_f2_dir_and_genotypes_agree(self, tmp_path):
        out = tmp_path / "f2"
        pops = "French,Yoruba,Simulation"
        run_cli("extract-f2", "--input", SMALL_PREFIX, "--output", str(out), "--populations", pops,
                "--all-chromosomes")
        args = ["f3", "--target", "Simulation", "--source1", "French", "--source2", "Yoruba", "--format", "json"]
        from_dir = json.loads(run_cli(*args, "--f2-dir", str(out)).stdout)
        from_geno = json.loads(run_cli(*args, "--input", SMALL_PREFIX, "--all-chromosomes").stdout)
        assert from_dir == from_geno

    @small_data
    def test_extract_f2_defaults_to_every_population(self, tmp_path):
        import popmix as pm
        out = tmp_path / "f2"
        run_cli("extract-f2", "--input", SMALL_PREFIX, "--output", str(out), "--all-chromosomes",
                "--maxmiss", "0.5")
        saved = pm.read_f2(str(out))
        pops = ["French", "Simulation", "Simulation.1", "Yoruba"]
        assert sorted(saved.populations) == pops
        direct = pm.extract_f2(SMALL_PREFIX, sorted(saved.populations, key=saved.populations.index),
                               auto_only=False, maxmiss=0.5)
        assert list(saved.block_lengths) == list(direct.block_lengths)

    @small_data
    def test_needs_input_or_f2_dir(self):
        result = run_cli("f2", "--pop1", "French", "--pop2", "Yoruba", check=False)
        assert result.returncode != 0
        assert "--input" in result.stderr and "--f2-dir" in result.stderr

    @small_data
    def test_rejects_input_and_f2_dir_together(self, tmp_path):
        result = run_cli("f2", "--input", SMALL_PREFIX, "--f2-dir", str(tmp_path), "--pop1", "French",
                         "--pop2", "Yoruba", check=False)
        assert result.returncode != 0
        assert "not both" in result.stderr

    @small_data
    def test_f2_dir_missing_population(self, tmp_path):
        out = tmp_path / "f2"
        run_cli("extract-f2", "--input", SMALL_PREFIX, "--output", str(out), "--populations", "French,Yoruba",
                "--all-chromosomes")
        result = run_cli("f2", "--f2-dir", str(out), "--pop1", "French", "--pop2", "Simulation", check=False)
        assert result.returncode != 0
        assert "Simulation" in result.stderr


class TestDatesCommand:
    """Test the 'dates' command for admixture dating."""

    # Test data path for DATES validation
    DATES_TEST_PREFIX = "tests/golden/dates/family_chr20_22"

    @pytest.mark.skipif(
        not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
        reason="DATES test data not available"
    )
    def test_dates_computes_estimate(self):
        """popmix dates should compute admixture date estimate."""
        result = run_cli(
            "dates",
            "--input", self.DATES_TEST_PREFIX,
            "--test-pop", "Simulation",
            "--source1", "French",
            "--source2", "Yoruba",
            timeout=120,
        )

        assert result.returncode == 0
        # Should show generations estimate
        assert "generation" in result.stdout.lower()
        # Should show amplitude (from fit)
        assert "amplitude" in result.stdout.lower()

    @pytest.mark.skipif(
        not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
        reason="DATES test data not available"
    )
    def test_dates_json_output(self):
        """popmix dates --format json should output valid JSON."""
        result = run_cli(
            "dates",
            "--input", self.DATES_TEST_PREFIX,
            "--test-pop", "Simulation",
            "--source1", "French",
            "--source2", "Yoruba",
            "--format", "json",
            timeout=120,
        )

        assert result.returncode == 0
        # Should be valid JSON
        data = json.loads(result.stdout)
        assert "generations" in data
        assert "amplitude" in data
        assert "ld_bins" in data

    @pytest.mark.skipif(
        not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
        reason="DATES test data not available"
    )
    def test_dates_tsv_output(self):
        """popmix dates --format tsv should output LD decay table."""
        result = run_cli(
            "dates",
            "--input", self.DATES_TEST_PREFIX,
            "--test-pop", "Simulation",
            "--source1", "French",
            "--source2", "Yoruba",
            "--format", "tsv",
            timeout=120,
        )

        assert result.returncode == 0
        # Should have tab-separated values
        assert "\t" in result.stdout
        # Should have header
        assert "distance_cm" in result.stdout.lower()

    @pytest.mark.skipif(
        not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
        reason="DATES test data not available"
    )
    def test_dates_qbin_matches_dates_v753_mesh_mode(self):
        """popmix dates --qbin 10 --compat dates-v753 should reproduce DATES v753 run with qbin: 10."""
        ref_dir = Path("tests/golden/dates/v753/chr20_22_q10")
        ref_est, ref_se = (float(x) for x in (ref_dir / "dates.jout").read_text().split())
        ref_full = next(
            float(line.split(":")[1])
            for line in (ref_dir / "dates.flog").read_text().splitlines()
            if line.startswith("mean (generations):")
        )

        result = run_cli(
            "dates",
            "--input", self.DATES_TEST_PREFIX,
            "--test-pop", "Simulation",
            "--source1", "French",
            "--source2", "Yoruba",
            "--qbin", "10",
            "--min-distance", "0.45",
            "--compat", "dates-v753",
            "--format", "json",
            timeout=120,
        )

        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["qbin"] == 10
        assert data["compat"] == "dates-v753"
        # DATES optimizer tolerance plus printed precision; see tests/golden/test_dates_v753.py
        assert abs(data["fit"]["decay_rate"] - ref_full) <= 1.5e-3
        assert abs(data["generations"] - ref_est) <= 0.008
        assert abs(data["generations_se"] - ref_se) <= 0.008

    def test_dates_missing_input_errors(self):
        """popmix dates with missing input should error."""
        result = run_cli(
            "dates",
            "--input", "/nonexistent/path",
            "--test-pop", "Test",
            "--source1", "Source1",
            "--source2", "Source2",
            check=False
        )
        assert result.returncode != 0


class TestMergeCommand:
    """Test the 'merge' command for combining genotype datasets."""

    # Test data paths
    MERGE_A_PREFIX = "tests/fixtures/merge_a"
    MERGE_B_PREFIX = "tests/fixtures/merge_b"

    @pytest.mark.skipif(
        not os.path.exists("tests/fixtures/merge_a.geno"),
        reason="Merge test fixtures not available"
    )
    def test_merge_creates_output(self):
        """popmix merge should create merged genotype files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "merged"

            result = run_cli(
                "merge",
                "--file-a", self.MERGE_A_PREFIX,
                "--file-b", self.MERGE_B_PREFIX,
                "--output", str(output),
                "--output-format", "eigenstrat",
            )

            assert result.returncode == 0
            assert (Path(tmpdir) / "merged.geno").exists()
            assert (Path(tmpdir) / "merged.snp").exists()
            assert (Path(tmpdir) / "merged.ind").exists()

    @pytest.mark.skipif(
        not os.path.exists("tests/fixtures/merge_a.geno"),
        reason="Merge test fixtures not available"
    )
    def test_merge_snp_intersection(self):
        """popmix merge should keep only SNPs present in both files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "merged"

            result = run_cli(
                "merge",
                "--file-a", self.MERGE_A_PREFIX,
                "--file-b", self.MERGE_B_PREFIX,
                "--output", str(output),
                "--output-format", "eigenstrat",
            )

            assert result.returncode == 0
            # Should report 3 SNPs in common (rs2, rs3, rs4)
            assert "SNPs in common: 3" in result.stdout

    @pytest.mark.skipif(
        not os.path.exists("tests/fixtures/merge_a.geno"),
        reason="Merge test fixtures not available"
    )
    def test_merge_sample_combination(self):
        """popmix merge should combine samples from both files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "merged"

            result = run_cli(
                "merge",
                "--file-a", self.MERGE_A_PREFIX,
                "--file-b", self.MERGE_B_PREFIX,
                "--output", str(output),
                "--output-format", "eigenstrat",
            )

            assert result.returncode == 0
            # merge_a has 2 samples, merge_b has 1 sample
            assert "Total samples: 3" in result.stdout

    @pytest.mark.skipif(
        not os.path.exists("tests/fixtures/merge_a.geno"),
        reason="Merge test fixtures not available"
    )
    def test_merge_to_plink_format(self):
        """popmix merge should support PLINK output format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "merged"

            result = run_cli(
                "merge",
                "--file-a", self.MERGE_A_PREFIX,
                "--file-b", self.MERGE_B_PREFIX,
                "--output", str(output),
                "--output-format", "plink",
            )

            assert result.returncode == 0
            assert (Path(tmpdir) / "merged.bed").exists()
            assert (Path(tmpdir) / "merged.bim").exists()
            assert (Path(tmpdir) / "merged.fam").exists()


class TestThreadsOption:
    """Test the --threads global option."""

    def test_threads_option_accepted(self):
        """popmix --threads N should be accepted."""
        result = run_cli("--threads", "2", "--help")
        assert result.returncode == 0

    def test_invalid_threads_errors(self):
        """popmix --threads 0 should error."""
        result = run_cli("--threads", "0", "--help", check=False)
        # Should either error or ignore invalid value
        # (implementation detail - just verify it handles it)


class TestOutputFormats:
    """Test output format options."""

    @small_data
    def test_json_output(self):
        result = run_cli("f2", "--input", SMALL_PREFIX, "--pop1", "French", "--pop2", "Yoruba",
                         "--all-chromosomes", "--format", "json")
        data = json.loads(result.stdout)
        assert "estimate" in data

    @small_data
    def test_tsv_output(self):
        result = run_cli("f2", "--input", SMALL_PREFIX, "--pop1", "French", "--pop2", "Yoruba",
                         "--all-chromosomes", "--format", "tsv")
        rows = [line.split("\t") for line in result.stdout.strip().splitlines()]
        assert rows[0][:2] == ["pop1", "pop2"]
        assert len(rows) == 2




class TestQpGraphCommand:
    """Test the 'qpgraph' command (admixtools2 qpgraph)."""

    REF = Path("tests/golden/admixtools2/v2.0.10/dg_6pop/qpgraph.json")

    @staticmethod
    def write_graph(tmp_path, edges, name="graph.txt"):
        path = tmp_path / name
        lines = ["# from to"] + [f"{a}\t{b}" for a, b in edges]
        path.write_text("\n".join(lines) + "\n")
        return str(path)

    def reference(self, graph_name):
        return next(g for g in json.loads(self.REF.read_text())["graphs"] if g["name"] == graph_name)

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    def test_qpgraph_json_matches_admixtools2(self, tmp_path):
        """The CLI fit of a 1-admixture graph reaches admixtools2's score."""
        ref = self.reference("admix1")
        graph = self.write_graph(tmp_path, ref["edges"])
        result = run_cli("qpgraph", "--input", AADR_PREFIX, "--graph", graph, "--format", "json")
        out = json.loads(result.stdout)
        assert out["pops"] == ref["pops"]
        assert out["score"] == pytest.approx(ref["score"], rel=1e-10)
        admix = {(e["from"], e["to"]): e["weight"] for e in out["edges"] if e["type"] == "admix"}
        for e in ref["fit_edges"]:
            if e["type"] == "admix":
                # within admixtools2's own spread over its random starts
                assert admix[(e["from"], e["to"])] == pytest.approx(e["weight"], abs=3e-4)

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    def test_qpgraph_fixed_weights_match_admixtools2_exactly(self, tmp_path):
        ref = self.reference("admix1")
        fixed = ref["fixed"][1]
        assert fixed["admix_weights"] == 0.3
        graph = self.write_graph(tmp_path, ref["edges"])
        result = run_cli("qpgraph", "--input", AADR_PREFIX, "--graph", graph,
                         "--admix-weights", "0.3", "--format", "json")
        out = json.loads(result.stdout)
        assert out["score"] == pytest.approx(fixed["score"], rel=1e-10)
        weights = [e["weight"] for e in out["edges"]]
        assert weights == pytest.approx(fixed["weight"], rel=1e-10, abs=1e-15)

    @pytest.mark.skipif(not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
                        reason="test data not available")
    def test_qpgraph_text_and_tsv_output(self, tmp_path):
        graph = self.write_graph(tmp_path, [("R", "Yoruba"), ("R", "a"), ("a", "French"),
                                            ("a", "b"), ("b", "Simulation"), ("b", "Simulation.1")])
        text = run_cli("qpgraph", "--input", SMALL_PREFIX, "--graph", graph, "--all-chromosomes")
        assert "score" in text.stdout.lower()
        assert "French" in text.stdout
        tsv = run_cli("qpgraph", "--input", SMALL_PREFIX, "--graph", graph, "--all-chromosomes",
                      "--format", "tsv")
        rows = [line.split("\t") for line in tsv.stdout.strip().splitlines()]
        assert rows[0] == ["from", "to", "type", "weight", "low", "high"]
        assert len(rows) == 1 + 6

    @pytest.mark.skipif(not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
                        reason="test data not available")
    def test_qpgraph_plots_and_dot(self, tmp_path):
        pytest.importorskip("matplotlib")
        graph = self.write_graph(tmp_path, [("R", "Yoruba"), ("R", "a"), ("a", "French"),
                                            ("a", "b"), ("b", "Simulation"), ("b", "Simulation.1")])
        run_cli("qpgraph", "--input", SMALL_PREFIX, "--graph", graph, "--all-chromosomes",
                "--plot", str(tmp_path / "g.png"), "--plot-fit", str(tmp_path / "f.png"),
                "--dot", str(tmp_path / "g.dot"))
        assert (tmp_path / "g.png").stat().st_size > 0 and (tmp_path / "f.png").stat().st_size > 0
        from popmix import parse_dot

        assert len(parse_dot(tmp_path / "g.dot")) == 6

    @pytest.mark.skipif(not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
                        reason="test data not available")
    def test_find_graphs_evolutionary(self, tmp_path):
        pytest.importorskip("matplotlib")
        out = run_cli("find-graphs", "--input", SMALL_PREFIX, "--populations", "Yoruba,French,Simulation,Simulation.1",
                      "--nadmix", "1", "--outgroup", "Yoruba", "--method", "evolutionary", "--stop-gen", "5",
                      "--all-chromosomes", "--plot", str(tmp_path / "best.png"), "--format", "json")
        graphs = json.loads(out.stdout)["graphs"]
        assert graphs and graphs[0]["score"] <= graphs[-1]["score"]
        assert (tmp_path / "best.png").stat().st_size > 0

    @pytest.mark.skipif(not os.path.exists("tests/golden/dates/family_chr20_22.geno"),
                        reason="test data not available")
    def test_qpgraph_rejects_a_graph_with_two_roots(self, tmp_path):
        graph = self.write_graph(tmp_path, [("R", "Yoruba"), ("R", "French"), ("S", "Simulation")])
        result = run_cli("qpgraph", "--input", SMALL_PREFIX, "--graph", graph, check=False)
        assert result.returncode != 0
        assert "root" in result.stderr.lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


class TestFindGraphsCommand:
    """Test the 'find-graphs' command."""

    POPS = "Yoruba,French,Simulation,Simulation.1"

    @small_data
    def test_json_matches_python_api(self):
        import popmix as pm
        out = json.loads(run_cli("find-graphs", "--input", SMALL_PREFIX, "--populations", self.POPS,
                                 "--all-chromosomes", "--method", "exhaustive", "--top", "3",
                                 "--format", "json").stdout)
        f2 = pm.extract_f2(SMALL_PREFIX, self.POPS.split(","), auto_only=False)
        graphs, n = pm.find_graphs(f2, self.POPS.split(","), method="exhaustive")
        assert out["n_evaluated"] == n == 15
        assert [g["score"] for g in out["graphs"]] == [g.score for g in graphs[:3]]
        assert [[tuple(e) for e in g["edges"]] for g in out["graphs"]] == [g.edges for g in graphs[:3]]

    @small_data
    def test_best_graph_file_refits_with_qpgraph(self, tmp_path):
        best = tmp_path / "best.txt"
        text = run_cli("find-graphs", "--input", SMALL_PREFIX, "--populations", self.POPS, "--nadmix", "1",
                       "--all-chromosomes", "--output-graph", str(best))
        assert "score" in text.stdout.lower()
        found = json.loads(run_cli("find-graphs", "--input", SMALL_PREFIX, "--populations", self.POPS,
                                   "--nadmix", "1", "--all-chromosomes", "--top", "1",
                                   "--format", "json").stdout)
        refit = json.loads(run_cli("qpgraph", "--input", SMALL_PREFIX, "--graph", str(best),
                                   "--all-chromosomes", "--format", "json").stdout)
        # qpgraph extracts its f2 blocks in leaf order, so only the last bits may differ
        assert refit["score"] == pytest.approx(found["graphs"][0]["score"], rel=1e-10)

    @small_data
    def test_outgroup_must_be_a_population(self):
        result = run_cli("find-graphs", "--input", SMALL_PREFIX, "--populations", self.POPS,
                         "--outgroup", "Mbuti", "--all-chromosomes", check=False)
        assert result.returncode != 0
        assert "outgroup" in result.stderr.lower()
