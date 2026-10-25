"""Tests for MultiQC custom module.

These tests verify the MultiQC integration for popmix result aggregation.
Following TDD: these tests are written BEFORE the implementation.
"""

import json

# =============================================================================
# Tests for File Pattern Matching
# =============================================================================

class TestFilePatternMatching:
    """Tests for MultiQC file pattern detection."""

    def test_detect_popmix_json_files(self, tmp_path):
        """Should detect *_popmix_stats.json files."""
        from popmix.multiqc.patterns import find_popmix_files

        # Create test files
        (tmp_path / "sample1_popmix_stats.json").write_text('{"f2": {}}')
        (tmp_path / "sample2_popmix_stats.json").write_text('{"f2": {}}')
        (tmp_path / "other.json").write_text("{}")

        files = find_popmix_files(str(tmp_path))

        assert len(files) == 2
        assert any("sample1" in f for f in files)
        assert any("sample2" in f for f in files)

    def test_detect_qpadm_json_files(self, tmp_path):
        """Should detect *_qpadm.json files."""
        from popmix.multiqc.patterns import find_popmix_files

        (tmp_path / "french_qpadm.json").write_text('{"target": "French"}')

        files = find_popmix_files(str(tmp_path))

        assert len(files) == 1
        assert "qpadm" in files[0]

    def test_detect_dates_json_files(self, tmp_path):
        """Should detect *_dates.json files."""
        from popmix.multiqc.patterns import find_popmix_files

        (tmp_path / "french_dates.json").write_text('{"generations": 150}')

        files = find_popmix_files(str(tmp_path))

        assert len(files) == 1

    def test_ignore_non_popmix_files(self, tmp_path):
        """Should ignore non-popmix JSON files."""
        from popmix.multiqc.patterns import find_popmix_files

        (tmp_path / "config.json").write_text("{}")
        (tmp_path / "data.json").write_text("{}")

        files = find_popmix_files(str(tmp_path))

        assert len(files) == 0


# =============================================================================
# Tests for File Parsing
# =============================================================================

class TestFileParsing:
    """Tests for parsing popmix output files."""

    def test_parse_f2_stats_json(self, tmp_path):
        """Should parse F2 statistics JSON."""
        from popmix.multiqc.parsers import parse_popmix_json

        data = {
            "type": "f2",
            "pop1": "French",
            "pop2": "Yoruba",
            "estimate": 0.05,
            "std_err": 0.002,
            "z_score": 25.0,
        }

        test_file = tmp_path / "test_popmix_stats.json"
        test_file.write_text(json.dumps(data))

        result = parse_popmix_json(str(test_file))

        assert result["type"] == "f2"
        assert result["estimate"] == 0.05

    def test_parse_qpadm_json(self, tmp_path):
        """Should parse qpAdm JSON."""
        from popmix.multiqc.parsers import parse_popmix_json

        data = {
            "target": "French",
            "sources": ["Yamnaya", "Anatolia_N"],
            "weights": {"Yamnaya": 0.45, "Anatolia_N": 0.55},
            "p_value": 0.237,
        }

        test_file = tmp_path / "test_qpadm.json"
        test_file.write_text(json.dumps(data))

        result = parse_popmix_json(str(test_file))

        assert result["target"] == "French"
        assert result["p_value"] == 0.237

    def test_parse_dates_json(self, tmp_path):
        """Should parse DATES JSON."""
        from popmix.multiqc.parsers import parse_popmix_json

        data = {
            "test_pop": "French",
            "generations": 150.2,
            "years_bp": 4356.0,
            "r_squared": 0.92,
        }

        test_file = tmp_path / "test_dates.json"
        test_file.write_text(json.dumps(data))

        result = parse_popmix_json(str(test_file))

        assert result["generations"] == 150.2


# =============================================================================
# Tests for Data Aggregation
# =============================================================================

class TestDataAggregation:
    """Tests for aggregating multiple sample results."""

    def test_aggregate_f2_results(self):
        """Should aggregate F2 results from multiple samples."""
        from popmix.multiqc.aggregator import aggregate_results

        results = [
            {"sample": "A", "type": "f2", "estimate": 0.05},
            {"sample": "B", "type": "f2", "estimate": 0.06},
        ]

        aggregated = aggregate_results(results)

        assert "f2" in aggregated
        assert len(aggregated["f2"]) == 2

    def test_aggregate_qpadm_results(self):
        """Should aggregate qpAdm results."""
        from popmix.multiqc.aggregator import aggregate_results

        results = [
            {"sample": "French", "type": "qpadm", "p_value": 0.24},
            {"sample": "Spanish", "type": "qpadm", "p_value": 0.15},
        ]

        aggregated = aggregate_results(results)

        assert "qpadm" in aggregated


# =============================================================================
# Tests for General Statistics
# =============================================================================

class TestGeneralStatistics:
    """Tests for MultiQC General Statistics columns."""

    def test_generate_general_stats_columns(self):
        """Should generate General Statistics column definitions."""
        from popmix.multiqc.general_stats import get_general_stats_headers

        headers = get_general_stats_headers()

        assert isinstance(headers, dict)
        assert len(headers) > 0

    def test_general_stats_include_pvalue(self):
        """General stats should include p-value column."""
        from popmix.multiqc.general_stats import get_general_stats_headers

        headers = get_general_stats_headers()

        # Should have a p-value related column
        assert any("p" in key.lower() for key in headers)

    def test_format_general_stats_data(self):
        """Should format data for General Statistics table."""
        from popmix.multiqc.general_stats import format_general_stats_data

        results = {
            "sample1": {"p_value": 0.24, "generations": 150},
        }

        data = format_general_stats_data(results)

        assert "sample1" in data
        assert isinstance(data["sample1"], dict)


# =============================================================================
# Tests for Search Patterns
# =============================================================================

class TestSearchPatterns:
    """Tests for MultiQC search patterns configuration."""

    def test_get_search_patterns(self):
        """Should return search patterns dictionary."""
        from popmix.multiqc.patterns import get_search_patterns

        patterns = get_search_patterns()

        assert isinstance(patterns, dict)
        assert "popmix" in patterns

    def test_search_patterns_have_required_keys(self):
        """Search patterns should have fn and contents keys."""
        from popmix.multiqc.patterns import get_search_patterns

        patterns = get_search_patterns()

        for pattern in patterns.values():
            assert "fn" in pattern or "contents" in pattern


# =============================================================================
# Tests for Module Integration
# =============================================================================

class TestModuleIntegration:
    """Tests for the main MultiQC module class."""

    def test_module_has_name(self):
        """Module should have a name attribute."""
        from popmix.multiqc import PopmixModule

        assert hasattr(PopmixModule, 'name') or hasattr(PopmixModule, 'get_name')

    def test_module_has_anchor(self):
        """Module should have an anchor for navigation."""
        from popmix.multiqc import PopmixModule

        assert hasattr(PopmixModule, 'anchor') or hasattr(PopmixModule, 'get_anchor')
