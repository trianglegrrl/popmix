"""Tests for batch analysis mode.

These tests verify the batch analysis commands for multi-target qpAdm and DATES.
Following TDD: these tests are written BEFORE the implementation.
"""

import json

# =============================================================================
# Tests for Batch Configuration Parsing
# =============================================================================

class TestBatchConfigParsing:
    """Tests for batch configuration file parsing."""

    def test_parse_targets_file(self, tmp_path):
        """Should parse targets from text file."""
        from popmix.cli.batch import parse_populations_file

        targets_file = tmp_path / "targets.txt"
        targets_file.write_text("French\nSpanish\nItalian\n")

        targets = parse_populations_file(str(targets_file))

        assert len(targets) == 3
        assert "French" in targets
        assert "Spanish" in targets
        assert "Italian" in targets

    def test_parse_populations_ignores_comments(self, tmp_path):
        """Should ignore comment lines starting with #."""
        from popmix.cli.batch import parse_populations_file

        file = tmp_path / "pops.txt"
        file.write_text("# This is a comment\nFrench\n# Another comment\nSpanish\n")

        pops = parse_populations_file(str(file))

        assert len(pops) == 2
        assert "French" in pops
        assert "Spanish" in pops

    def test_parse_populations_ignores_empty_lines(self, tmp_path):
        """Should ignore empty lines."""
        from popmix.cli.batch import parse_populations_file

        file = tmp_path / "pops.txt"
        file.write_text("French\n\nSpanish\n\n\nItalian\n")

        pops = parse_populations_file(str(file))

        assert len(pops) == 3

    def test_parse_yaml_config(self, tmp_path):
        """Should parse YAML configuration file."""
        from popmix.cli.batch import parse_batch_config

        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
targets:
  - French
  - Spanish
sources:
  - Yamnaya
  - Anatolia_N
  - WHG
outgroups:
  - Mbuti
  - Han
""")

        config = parse_batch_config(str(config_file))

        assert len(config["targets"]) == 2
        assert len(config["sources"]) == 3
        assert len(config["outgroups"]) == 2


# =============================================================================
# Tests for Batch Job Generation
# =============================================================================

class TestBatchJobGeneration:
    """Tests for generating batch analysis jobs."""

    def test_generate_qpadm_jobs(self):
        """Should generate qpAdm jobs for all targets."""
        from popmix.cli.batch import generate_qpadm_jobs

        targets = ["French", "Spanish"]
        sources = ["Yamnaya", "Anatolia_N"]
        outgroups = ["Mbuti", "Han"]

        jobs = generate_qpadm_jobs(targets, sources, outgroups)

        assert len(jobs) == 2  # One per target
        assert jobs[0]["target"] == "French"
        assert jobs[1]["target"] == "Spanish"

    def test_qpadm_job_includes_sources_and_outgroups(self):
        """Each qpAdm job should include sources and outgroups."""
        from popmix.cli.batch import generate_qpadm_jobs

        jobs = generate_qpadm_jobs(
            targets=["French"],
            sources=["Yamnaya", "Anatolia_N"],
            outgroups=["Mbuti"],
        )

        job = jobs[0]
        assert "sources" in job
        assert "outgroups" in job
        assert len(job["sources"]) == 2

    def test_generate_dates_jobs(self):
        """Should generate DATES jobs for all targets."""
        from popmix.cli.batch import generate_dates_jobs

        targets = ["French", "Spanish"]
        source1 = "Yamnaya"
        source2 = "Anatolia_N"

        jobs = generate_dates_jobs(targets, source1, source2)

        assert len(jobs) == 2
        assert jobs[0]["test_pop"] == "French"
        assert jobs[0]["source1"] == "Yamnaya"
        assert jobs[0]["source2"] == "Anatolia_N"


# =============================================================================
# Tests for Batch Execution
# =============================================================================

class TestBatchExecution:
    """Tests for batch job execution."""

    def test_batch_runner_tracks_progress(self):
        """Batch runner should track completed jobs."""
        from popmix.cli.batch import BatchRunner

        runner = BatchRunner(n_jobs=3)

        assert runner.total == 3
        assert runner.completed == 0

    def test_batch_runner_updates_on_completion(self):
        """Runner should update completed count."""
        from popmix.cli.batch import BatchRunner

        runner = BatchRunner(n_jobs=3)
        runner.mark_completed()

        assert runner.completed == 1

    def test_batch_runner_calculates_progress(self):
        """Runner should calculate progress percentage."""
        from popmix.cli.batch import BatchRunner

        runner = BatchRunner(n_jobs=4)
        runner.mark_completed()
        runner.mark_completed()

        assert runner.progress_percent == 50.0


# =============================================================================
# Tests for Result Collection
# =============================================================================

class TestResultCollection:
    """Tests for collecting batch results."""

    def test_collect_results_returns_list(self):
        """Result collector should return list of results."""
        from popmix.cli.batch import ResultCollector

        collector = ResultCollector()
        collector.add_result("French", {"p_value": 0.24})
        collector.add_result("Spanish", {"p_value": 0.15})

        results = collector.get_results()

        assert len(results) == 2

    def test_collect_results_includes_metadata(self):
        """Each result should include target name."""
        from popmix.cli.batch import ResultCollector

        collector = ResultCollector()
        collector.add_result("French", {"p_value": 0.24})

        results = collector.get_results()

        assert results[0]["target"] == "French"

    def test_collect_results_tracks_failures(self):
        """Should track failed analyses."""
        from popmix.cli.batch import ResultCollector

        collector = ResultCollector()
        collector.add_result("French", {"p_value": 0.24})
        collector.add_failure("Spanish", "Insufficient outgroups")

        assert collector.n_success == 1
        assert collector.n_failures == 1


# =============================================================================
# Tests for Output Generation
# =============================================================================

class TestBatchOutputGeneration:
    """Tests for batch output file generation."""

    def test_write_batch_results_json(self, tmp_path):
        """Should write all results to single JSON file."""
        from popmix.cli.batch import write_batch_results

        results = [
            {"target": "French", "p_value": 0.24},
            {"target": "Spanish", "p_value": 0.15},
        ]
        output_path = tmp_path / "batch_results.json"

        write_batch_results(results, str(output_path), format="json")

        assert output_path.exists()
        with open(output_path) as f:
            data = json.load(f)
        assert len(data) == 2

    def test_write_batch_results_creates_per_target_files(self, tmp_path):
        """Should create individual result files per target."""
        from popmix.cli.batch import write_batch_results_individual

        results = [
            {"target": "French", "p_value": 0.24},
            {"target": "Spanish", "p_value": 0.15},
        ]

        write_batch_results_individual(results, str(tmp_path), format="json")

        assert (tmp_path / "French_qpadm.json").exists()
        assert (tmp_path / "Spanish_qpadm.json").exists()


# =============================================================================
# Tests for Summary Report
# =============================================================================

class TestBatchSummaryReport:
    """Tests for batch summary report generation."""

    def test_generate_summary_statistics(self):
        """Should generate summary statistics."""
        from popmix.cli.batch import generate_batch_summary

        results = [
            {"target": "French", "p_value": 0.24, "feasible": True},
            {"target": "Spanish", "p_value": 0.03, "feasible": True},
            {"target": "Italian", "p_value": 0.001, "feasible": False},
        ]

        summary = generate_batch_summary(results)

        assert summary["n_total"] == 3
        assert summary["n_good_fit"] == 1  # p > 0.05
        assert summary["n_feasible"] == 2

    def test_summary_includes_best_models(self):
        """Summary should include best fitting models."""
        from popmix.cli.batch import generate_batch_summary

        results = [
            {"target": "French", "p_value": 0.50},
            {"target": "Spanish", "p_value": 0.03},
        ]

        summary = generate_batch_summary(results)

        assert "best_models" in summary
        assert summary["best_models"][0]["target"] == "French"  # Highest p-value
