"""Tests for output formatters.

These tests verify the machine-readable output formats for popmix results.
Following TDD: these tests are written BEFORE the implementation.
"""

import json

# =============================================================================
# Mock result classes (reused from test_rich_output.py)
# =============================================================================

class MockF2Result:
    """Mock F2 result for testing."""

    def __init__(
        self,
        pop1: str = "PopA",
        pop2: str = "PopB",
        estimate: float = 0.02,
        std_err: float = 0.001,
        z_score: float = 20.0,
        n_snps: int = 100000,
    ):
        self.pop1 = pop1
        self.pop2 = pop2
        self.estimate = estimate
        self.std_err = std_err
        self.z_score = z_score
        self.n_snps = n_snps


class MockQpAdmResult:
    """Mock qpAdm result for testing."""

    def __init__(
        self,
        target: str = "French",
        sources: list = None,
        weights: list = None,
        std_errs: list = None,
        chi_sq: float = 4.23,
        df: int = 3,
        p_value: float = 0.237,
    ):
        self.target = target
        self.sources = sources or ["Yamnaya", "Anatolia_N"]
        self.weights = weights or [0.45, 0.55]
        self.std_errs = std_errs or [0.036, 0.036]
        self.chi_sq = chi_sq
        self.df = df
        self.p_value = p_value

    def is_feasible(self) -> bool:
        return all(0 <= w <= 1 for w in self.weights)

    def is_good_fit(self) -> bool:
        return self.p_value > 0.05


class MockDatesResult:
    """Mock DATES result for testing."""

    def __init__(
        self,
        test_pop: str = "French",
        source1: str = "Yamnaya",
        source2: str = "Anatolia_N",
        generations: float = 150.2,
        generations_se: float = 12.3,
        years_bp: float = 4356.0,
        years_bp_se: float = 357.0,
        amplitude: float = 0.036,
        r_squared: float = 0.92,
        n_snps: int = 1200000,
    ):
        self.test_pop = test_pop
        self.source1 = source1
        self.source2 = source2
        self.generations = generations
        self.generations_se = generations_se
        self.years_bp = years_bp
        self.years_bp_se = years_bp_se
        self.amplitude = amplitude
        self.r_squared = r_squared
        self.n_snps = n_snps


# =============================================================================
# Tests for JSON formatting
# =============================================================================

class TestJsonFormatter:
    """Tests for JSON output formatting."""

    def test_format_f2_as_json_valid(self):
        """JSON output should be valid JSON."""
        from popmix.cli.formatters import format_as_json

        result = MockF2Result()
        output = format_as_json(result, "f2")

        # Should be valid JSON
        parsed = json.loads(output)
        assert isinstance(parsed, dict)

    def test_format_f2_as_json_contains_required_fields(self):
        """F2 JSON should contain all required fields."""
        from popmix.cli.formatters import format_as_json

        result = MockF2Result(
            pop1="A", pop2="B", estimate=0.02, std_err=0.001, z_score=20.0
        )
        output = format_as_json(result, "f2")
        parsed = json.loads(output)

        assert "pop1" in parsed
        assert "pop2" in parsed
        assert "estimate" in parsed
        assert "std_err" in parsed
        assert "z_score" in parsed

    def test_format_f2_as_json_includes_confidence_interval(self):
        """F2 JSON should include confidence interval."""
        from popmix.cli.formatters import format_as_json

        result = MockF2Result(estimate=0.02, std_err=0.001)
        output = format_as_json(result, "f2")
        parsed = json.loads(output)

        assert "ci_lower" in parsed or "confidence_interval" in parsed
        assert "ci_upper" in parsed or "confidence_interval" in parsed

    def test_format_qpadm_as_json_contains_weights(self):
        """qpAdm JSON should contain weights dictionary."""
        from popmix.cli.formatters import format_as_json

        result = MockQpAdmResult()
        output = format_as_json(result, "qpadm")
        parsed = json.loads(output)

        assert "weights" in parsed
        assert "target" in parsed
        assert "sources" in parsed

    def test_format_dates_as_json_contains_timing(self):
        """DATES JSON should contain timing estimates."""
        from popmix.cli.formatters import format_as_json

        result = MockDatesResult()
        output = format_as_json(result, "dates")
        parsed = json.loads(output)

        assert "generations" in parsed
        assert "years_bp" in parsed


# =============================================================================
# Tests for YAML formatting
# =============================================================================

class TestYamlFormatter:
    """Tests for YAML output formatting."""

    def test_format_f2_as_yaml_valid(self):
        """YAML output should be valid YAML."""
        import yaml

        from popmix.cli.formatters import format_as_yaml

        result = MockF2Result()
        output = format_as_yaml(result, "f2")

        # Should be valid YAML
        parsed = yaml.safe_load(output)
        assert isinstance(parsed, dict)

    def test_format_f2_as_yaml_contains_required_fields(self):
        """F2 YAML should contain all required fields."""
        import yaml

        from popmix.cli.formatters import format_as_yaml

        result = MockF2Result()
        output = format_as_yaml(result, "f2")
        parsed = yaml.safe_load(output)

        assert "pop1" in parsed
        assert "pop2" in parsed
        assert "estimate" in parsed

    def test_yaml_is_human_readable(self):
        """YAML output should be human-readable (not minified)."""
        from popmix.cli.formatters import format_as_yaml

        result = MockF2Result()
        output = format_as_yaml(result, "f2")

        # Should have multiple lines
        assert output.count("\n") > 3


# =============================================================================
# Tests for TSV formatting
# =============================================================================

class TestTsvFormatter:
    """Tests for TSV output formatting."""

    def test_format_f2_as_tsv_has_header(self):
        """TSV output should have a header row."""
        from popmix.cli.formatters import format_as_tsv

        result = MockF2Result()
        output = format_as_tsv(result, "f2")
        lines = output.strip().split("\n")

        assert len(lines) >= 2  # Header + data
        assert "pop1" in lines[0].lower() or "estimate" in lines[0].lower()

    def test_format_f2_as_tsv_tab_separated(self):
        """TSV output should use tab separators."""
        from popmix.cli.formatters import format_as_tsv

        result = MockF2Result()
        output = format_as_tsv(result, "f2")
        lines = output.strip().split("\n")

        # Each line should have tabs
        assert "\t" in lines[0]
        assert "\t" in lines[1]

    def test_format_f2_as_tsv_correct_columns(self):
        """TSV should have correct column count."""
        from popmix.cli.formatters import format_as_tsv

        result = MockF2Result()
        output = format_as_tsv(result, "f2")
        lines = output.strip().split("\n")

        header_cols = len(lines[0].split("\t"))
        data_cols = len(lines[1].split("\t"))

        assert header_cols == data_cols


# =============================================================================
# Tests for NDJSON formatting
# =============================================================================

class TestNdjsonFormatter:
    """Tests for NDJSON (newline-delimited JSON) output formatting."""

    def test_format_as_ndjson_one_line_per_result(self):
        """NDJSON should have one JSON object per line."""
        from popmix.cli.formatters import format_as_ndjson

        results = [
            MockF2Result(pop1="A", pop2="B"),
            MockF2Result(pop1="C", pop2="D"),
        ]
        output = format_as_ndjson(results, "f2")
        lines = output.strip().split("\n")

        assert len(lines) == 2

        # Each line should be valid JSON
        for line in lines:
            parsed = json.loads(line)
            assert isinstance(parsed, dict)

    def test_ndjson_each_line_valid_json(self):
        """Each line in NDJSON should be valid JSON."""
        from popmix.cli.formatters import format_as_ndjson

        results = [MockF2Result(), MockF2Result()]
        output = format_as_ndjson(results, "f2")

        for line in output.strip().split("\n"):
            # Should not raise
            json.loads(line)


# =============================================================================
# Tests for Parquet formatting
# =============================================================================

class TestParquetFormatter:
    """Tests for Parquet output formatting."""

    def test_write_f2_as_parquet_creates_file(self, tmp_path):
        """Parquet writer should create a file."""
        from popmix.cli.formatters import write_as_parquet

        result = MockF2Result()
        output_path = tmp_path / "test.parquet"

        write_as_parquet(result, "f2", str(output_path))

        assert output_path.exists()

    def test_parquet_readable_by_pandas(self, tmp_path):
        """Parquet file should be readable by pandas."""
        import pandas as pd

        from popmix.cli.formatters import write_as_parquet

        result = MockF2Result(pop1="A", pop2="B", estimate=0.02)
        output_path = tmp_path / "test.parquet"

        write_as_parquet(result, "f2", str(output_path))

        df = pd.read_parquet(output_path)
        assert len(df) >= 1

    def test_parquet_contains_correct_data(self, tmp_path):
        """Parquet file should contain correct data."""
        import pandas as pd

        from popmix.cli.formatters import write_as_parquet

        result = MockF2Result(pop1="TestA", pop2="TestB", estimate=0.0234)
        output_path = tmp_path / "test.parquet"

        write_as_parquet(result, "f2", str(output_path))

        df = pd.read_parquet(output_path)
        assert df.iloc[0]["pop1"] == "TestA"
        assert df.iloc[0]["pop2"] == "TestB"


# =============================================================================
# Tests for file output
# =============================================================================

class TestFileOutput:
    """Tests for writing results to files."""

    def test_write_to_json_file(self, tmp_path):
        """Should write JSON to file."""
        from popmix.cli.formatters import write_result_to_file

        result = MockF2Result()
        output_path = tmp_path / "result.json"

        write_result_to_file(result, "f2", str(output_path), format="json")

        assert output_path.exists()
        with open(output_path) as f:
            data = json.load(f)
        assert "estimate" in data

    def test_write_to_yaml_file(self, tmp_path):
        """Should write YAML to file."""
        import yaml

        from popmix.cli.formatters import write_result_to_file

        result = MockF2Result()
        output_path = tmp_path / "result.yaml"

        write_result_to_file(result, "f2", str(output_path), format="yaml")

        assert output_path.exists()
        with open(output_path) as f:
            data = yaml.safe_load(f)
        assert "estimate" in data

    def test_write_auto_detects_format_from_extension(self, tmp_path):
        """Should auto-detect format from file extension."""
        from popmix.cli.formatters import write_result_to_file

        result = MockF2Result()

        # JSON extension
        json_path = tmp_path / "result.json"
        write_result_to_file(result, "f2", str(json_path))
        assert json_path.exists()
        with open(json_path) as f:
            json.load(f)  # Should be valid JSON

        # YAML extension
        yaml_path = tmp_path / "result.yaml"
        write_result_to_file(result, "f2", str(yaml_path))
        assert yaml_path.exists()

    def test_write_creates_parent_directories(self, tmp_path):
        """Should create parent directories if they don't exist."""
        from popmix.cli.formatters import write_result_to_file

        result = MockF2Result()
        output_path = tmp_path / "nested" / "dir" / "result.json"

        write_result_to_file(result, "f2", str(output_path), format="json")

        assert output_path.exists()


# =============================================================================
# Tests for result serialization
# =============================================================================

class TestResultSerialization:
    """Tests for converting result objects to dictionaries."""

    def test_serialize_f2_result(self):
        """F2 result should serialize to dict."""
        from popmix.cli.formatters import result_to_dict

        result = MockF2Result(pop1="A", pop2="B", estimate=0.02)
        data = result_to_dict(result, "f2")

        assert isinstance(data, dict)
        assert data["pop1"] == "A"
        assert data["pop2"] == "B"
        assert data["estimate"] == 0.02

    def test_serialize_qpadm_result(self):
        """qpAdm result should serialize to dict with nested structure."""
        from popmix.cli.formatters import result_to_dict

        result = MockQpAdmResult()
        data = result_to_dict(result, "qpadm")

        assert "target" in data
        assert "weights" in data
        assert "model_fit" in data or "chi_sq" in data

    def test_serialize_includes_metadata(self):
        """Serialized results should include metadata."""
        from popmix.cli.formatters import result_to_dict

        result = MockF2Result()
        data = result_to_dict(result, "f2", include_metadata=True)

        assert "result_type" in data or "type" in data
