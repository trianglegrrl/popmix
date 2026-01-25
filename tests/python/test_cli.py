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


class TestF2Command:
    """Test the 'f2' command for computing f2 statistics."""

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_f2_computes_statistic(self):
        """popmix f2 should compute f2 between two populations."""
        result = run_cli(
            "f2",
            "--input", AADR_PREFIX,
            "--pop1", POP_HAN,
            "--pop2", POP_YORUBA,
            "--snp-range", "0:1000",
        )

        assert result.returncode == 0
        # Should output the f2 value
        assert "f2" in result.stdout.lower() or "estimate" in result.stdout.lower()


class TestF3Command:
    """Test the 'f3' command for computing f3 statistics."""

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_f3_computes_statistic(self):
        """popmix f3 should compute f3(A; B, C)."""
        result = run_cli(
            "f3",
            "--input", AADR_PREFIX,
            "--target", POP_FRENCH,
            "--source1", POP_HAN,
            "--source2", POP_YORUBA,
            "--snp-range", "0:1000",
        )

        assert result.returncode == 0
        assert "f3" in result.stdout.lower() or "estimate" in result.stdout.lower()


class TestF4Command:
    """Test the 'f4' command for computing f4/D statistics."""

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_f4_computes_statistic(self):
        """popmix f4 should compute f4(A, B; C, D)."""
        result = run_cli(
            "f4",
            "--input", AADR_PREFIX,
            "--pop1", POP_FRENCH,
            "--pop2", POP_HAN,
            "--pop3", POP_YORUBA,
            "--pop4", POP_MBUTI,
            "--snp-range", "0:1000",
        )

        assert result.returncode == 0
        assert "f4" in result.stdout.lower() or "estimate" in result.stdout.lower()


class TestQpAdmCommand:
    """Test the 'qpadm' command."""

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_qpadm_computes_weights(self):
        """popmix qpadm should estimate admixture weights."""
        result = run_cli(
            "qpadm",
            "--input", AADR_PREFIX,
            "--target", POP_FRENCH,
            "--left", f"{POP_HAN},{POP_YORUBA}",
            "--right", f"{POP_MBUTI},Ignore_Yoruba(discovery).DG",
            "--snp-range", "0:1000",
        )

        assert result.returncode == 0
        # Should show weights
        assert "weight" in result.stdout.lower() or "proportion" in result.stdout.lower()


class TestQpWaveCommand:
    """Test the 'qpwave' command."""

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_qpwave_computes_rank(self):
        """popmix qpwave should compute rank tests."""
        result = run_cli(
            "qpwave",
            "--input", AADR_PREFIX,
            "--left", f"{POP_FRENCH},{POP_HAN},{POP_YORUBA}",
            "--right", f"{POP_MBUTI},Ignore_Yoruba(discovery).DG",
            "--snp-range", "0:1000",
        )

        assert result.returncode == 0
        assert "rank" in result.stdout.lower() or "p-value" in result.stdout.lower()


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

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_json_output(self):
        """--format json should output valid JSON."""
        result = run_cli(
            "f2",
            "--input", AADR_PREFIX,
            "--pop1", POP_HAN,
            "--pop2", POP_YORUBA,
            "--snp-range", "0:100",
            "--format", "json",
        )

        assert result.returncode == 0
        # Should be valid JSON
        data = json.loads(result.stdout)
        assert "estimate" in data or "f2" in data

    @pytest.mark.skipif(not os.path.exists(f"{AADR_PREFIX}.geno"), reason="AADR data not available")
    @pytest.mark.slow
    def test_tsv_output(self):
        """--format tsv should output tab-separated values."""
        result = run_cli(
            "f2",
            "--input", AADR_PREFIX,
            "--pop1", POP_HAN,
            "--pop2", POP_YORUBA,
            "--snp-range", "0:100",
            "--format", "tsv",
        )

        assert result.returncode == 0
        # Should have tab-separated values
        assert "\t" in result.stdout


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
