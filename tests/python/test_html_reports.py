"""Tests for HTML report generation.

These tests verify the Jinja2-based HTML report generator with Plotly visualizations.
Following TDD: these tests are written BEFORE the implementation.
"""




# =============================================================================
# Mock result classes
# =============================================================================

class MockF2Result:
    def __init__(self):
        self.pop1 = "A"
        self.pop2 = "B"
        self.estimate = 0.02
        self.std_err = 0.001
        self.z_score = 20.0
        self.n_snps = 100000


class MockQpAdmResult:
    def __init__(self):
        self.target = "French"
        self.sources = ["Yamnaya", "Anatolia_N"]
        self.weights = [0.45, 0.55]
        self.std_errs = [0.036, 0.036]
        self.chi_sq = 4.23
        self.df = 3
        self.p_value = 0.237

    def is_feasible(self):
        return all(0 <= w <= 1 for w in self.weights)

    def is_good_fit(self):
        return self.p_value > 0.05


class MockDatesResult:
    def __init__(self):
        self.test_pop = "French"
        self.source1 = "Yamnaya"
        self.source2 = "Anatolia_N"
        self.generations = 150.2
        self.generations_se = 12.3
        self.years_bp = 4356.0
        self.years_bp_se = 357.0
        self.amplitude = 0.036
        self.r_squared = 0.92
        self.n_snps = 1200000
        self.ld_bins = [
            {"distance_cm": 0.1, "correlation": 0.137, "n_pairs": 298384},
            {"distance_cm": 0.2, "correlation": 0.083, "n_pairs": 206537},
            {"distance_cm": 0.3, "correlation": 0.050, "n_pairs": 150000},
        ]


# =============================================================================
# Tests for HTML Report Generation
# =============================================================================

class TestHtmlReportBasics:
    """Basic tests for HTML report generation."""

    def test_generate_html_report_returns_string(self):
        """Should return a string."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        assert isinstance(html, str)
        assert len(html) > 0

    def test_html_report_is_valid_html(self):
        """Generated HTML should be valid HTML5."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        assert "<!DOCTYPE html>" in html
        assert "<html" in html
        assert "</html>" in html
        assert "<head>" in html
        assert "</head>" in html
        assert "<body>" in html
        assert "</body>" in html

    def test_html_report_includes_title(self):
        """HTML report should include a title."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm", title="Test Report")

        assert "<title>" in html
        assert "Test Report" in html

    def test_html_report_includes_popmix_branding(self):
        """Report should include popmix branding."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        assert "popmix" in html.lower()


# =============================================================================
# Tests for Plotly Integration
# =============================================================================

class TestPlotlyIntegration:
    """Tests for Plotly.js chart integration."""

    def test_html_includes_plotly_js(self):
        """HTML should include Plotly.js library."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        assert "plotly" in html.lower()

    def test_qpadm_report_contains_bar_chart(self):
        """qpAdm report should contain stacked bar chart."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        # Should have chart div
        assert "chart" in html.lower() or "plot" in html.lower()

    def test_dates_report_contains_decay_curve(self):
        """DATES report should contain LD decay curve."""
        from popmix.reports.html_report import generate_html_report

        result = MockDatesResult()
        html = generate_html_report(result, "dates")

        # Should reference decay or LD
        assert "decay" in html.lower() or "ld" in html.lower()


# =============================================================================
# Tests for qpAdm Report Content
# =============================================================================

class TestQpAdmReportContent:
    """Tests for qpAdm-specific report content."""

    def test_qpadm_report_shows_target(self):
        """qpAdm report should show target population."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        result.target = "TestTarget"
        html = generate_html_report(result, "qpadm")

        assert "TestTarget" in html

    def test_qpadm_report_shows_sources(self):
        """qpAdm report should show source populations."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        result.sources = ["Source1", "Source2"]
        html = generate_html_report(result, "qpadm")

        assert "Source1" in html
        assert "Source2" in html

    def test_qpadm_report_shows_weights(self):
        """qpAdm report should show admixture weights."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        result.weights = [0.45, 0.55]
        html = generate_html_report(result, "qpadm")

        # Should show weights as percentages
        assert "45" in html or "0.45" in html
        assert "55" in html or "0.55" in html

    def test_qpadm_report_shows_model_fit(self):
        """qpAdm report should show model fit statistics."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        # Should show p-value and model fit indicator
        assert "p-value" in html.lower() or "0.237" in html
        assert "good fit" in html.lower() or "model fit" in html.lower()


# =============================================================================
# Tests for DATES Report Content
# =============================================================================

class TestDatesReportContent:
    """Tests for DATES-specific report content."""

    def test_dates_report_shows_timing(self):
        """DATES report should show timing estimates."""
        from popmix.reports.html_report import generate_html_report

        result = MockDatesResult()
        html = generate_html_report(result, "dates")

        assert "generation" in html.lower()
        assert "150" in html  # generations

    def test_dates_report_shows_years_bp(self):
        """DATES report should show years before present."""
        from popmix.reports.html_report import generate_html_report

        result = MockDatesResult()
        html = generate_html_report(result, "dates")

        assert "year" in html.lower() or "BP" in html

    def test_dates_report_shows_fit_quality(self):
        """DATES report should show R² fit quality."""
        from popmix.reports.html_report import generate_html_report

        result = MockDatesResult()
        html = generate_html_report(result, "dates")

        assert "R²" in html or "r-squared" in html.lower() or "0.92" in html


# =============================================================================
# Tests for File Output
# =============================================================================

class TestHtmlFileOutput:
    """Tests for writing HTML reports to files."""

    def test_write_report_creates_file(self, tmp_path):
        """Should create HTML file."""
        from popmix.reports.html_report import write_html_report

        result = MockQpAdmResult()
        output_path = tmp_path / "report.html"

        write_html_report(result, "qpadm", str(output_path))

        assert output_path.exists()

    def test_written_file_is_valid_html(self, tmp_path):
        """Written file should contain valid HTML."""
        from popmix.reports.html_report import write_html_report

        result = MockQpAdmResult()
        output_path = tmp_path / "report.html"

        write_html_report(result, "qpadm", str(output_path))

        content = output_path.read_text()
        assert "<!DOCTYPE html>" in content

    def test_write_creates_parent_directories(self, tmp_path):
        """Should create parent directories if needed."""
        from popmix.reports.html_report import write_html_report

        result = MockQpAdmResult()
        output_path = tmp_path / "nested" / "dir" / "report.html"

        write_html_report(result, "qpadm", str(output_path))

        assert output_path.exists()


# =============================================================================
# Tests for Styling
# =============================================================================

class TestHtmlStyling:
    """Tests for HTML report styling."""

    def test_html_includes_css(self):
        """HTML should include CSS styling."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        assert "<style>" in html or "stylesheet" in html.lower()

    def test_html_is_responsive(self):
        """HTML should include responsive meta tag."""
        from popmix.reports.html_report import generate_html_report

        result = MockQpAdmResult()
        html = generate_html_report(result, "qpadm")

        assert "viewport" in html
