"""Tests for Rich terminal output formatting.

These tests verify the Rich-based output formatting for popmix CLI commands.
Following TDD: these tests are written BEFORE the implementation.
"""



# Mock result classes for testing
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


class MockF3Result:
    """Mock F3 result for testing."""

    def __init__(
        self,
        target: str = "Target",
        source1: str = "Source1",
        source2: str = "Source2",
        estimate: float = -0.001,
        std_err: float = 0.0001,
        z_score: float = -10.0,
        n_blocks: int = 500,
    ):
        self.target = target
        self.source1 = source1
        self.source2 = source2
        self.estimate = estimate
        self.std_err = std_err
        self.z_score = z_score
        self.n_blocks = n_blocks


class MockF4Result:
    """Mock F4 result for testing."""

    def __init__(
        self,
        pop1: str = "A",
        pop2: str = "B",
        pop3: str = "C",
        pop4: str = "D",
        estimate: float = 0.0005,
        std_err: float = 0.0001,
        z_score: float = 5.0,
        n_blocks: int = 500,
    ):
        self.pop1 = pop1
        self.pop2 = pop2
        self.pop3 = pop3
        self.pop4 = pop4
        self.estimate = estimate
        self.std_err = std_err
        self.z_score = z_score
        self.n_blocks = n_blocks


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
        """Check if weights are in [0, 1]."""
        return all(0 <= w <= 1 for w in self.weights)

    def is_good_fit(self) -> bool:
        """Check if p-value > 0.05."""
        return self.p_value > 0.05


class MockQpWaveRankTest:
    """Mock qpWave rank test for testing."""

    def __init__(self, rank: int, chi_sq: float, df: int, p_value: float):
        self.rank = rank
        self.chi_sq = chi_sq
        self.df = df
        self.p_value = p_value


class MockQpWaveResult:
    """Mock qpWave result for testing."""

    def __init__(
        self,
        left_pops: list = None,
        right_pops: list = None,
        estimated_rank: int = 1,
        rank_tests: list = None,
    ):
        self.left_pops = left_pops or ["French", "Spanish", "Italian"]
        self.right_pops = right_pops or ["Mbuti", "Han", "Papuan"]
        self.estimated_rank = estimated_rank
        self.rank_tests = rank_tests or [
            MockQpWaveRankTest(0, 45.23, 6, 0.000001),
            MockQpWaveRankTest(1, 3.12, 2, 0.21),
            MockQpWaveRankTest(2, 0.0, 0, 1.0),
        ]

    def is_clade(self) -> bool:
        """Check if rank 0 is not rejected."""
        return self.rank_tests[0].p_value > 0.05


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
# Tests for F2 formatting
# =============================================================================

class TestF2Formatting:
    """Tests for F2 result formatting."""

    def test_format_f2_result_contains_estimate(self):
        """F2 output should contain the estimate value."""
        from popmix.cli.rich_output import format_f2_result

        result = MockF2Result(estimate=0.023456)
        output = format_f2_result(result)

        assert "0.0235" in output or "0.023456" in output

    def test_format_f2_result_contains_std_err(self):
        """F2 output should contain the standard error."""
        from popmix.cli.rich_output import format_f2_result

        result = MockF2Result(std_err=0.000123)
        output = format_f2_result(result)

        assert "0.000123" in output or "0.0001" in output

    def test_format_f2_result_contains_z_score(self):
        """F2 output should contain the Z-score."""
        from popmix.cli.rich_output import format_f2_result

        result = MockF2Result(z_score=190.5)
        output = format_f2_result(result)

        assert "190" in output

    def test_format_f2_result_contains_confidence_interval(self):
        """F2 output should contain 95% confidence interval."""
        from popmix.cli.rich_output import format_f2_result

        # 95% CI = estimate ± 1.96 * std_err
        # 0.02 ± 1.96 * 0.001 = [0.01804, 0.02196]
        result = MockF2Result(estimate=0.02, std_err=0.001)
        output = format_f2_result(result)

        assert "[" in output and "]" in output

    def test_format_f2_result_significant_colored_green(self):
        """Significant F2 result (|Z| >= 3) should indicate significance."""
        from popmix.cli.rich_output import format_f2_result

        result = MockF2Result(z_score=5.0)
        output = format_f2_result(result)

        # Should contain significance indicator
        assert "significant" in output.lower() or "✓" in output or "green" in output.lower()

    def test_format_f2_result_not_significant_gray(self):
        """Non-significant F2 result (|Z| < 2) should indicate non-significance."""
        from popmix.cli.rich_output import format_f2_result

        result = MockF2Result(z_score=1.5)
        output = format_f2_result(result)

        # Should not claim significance
        assert "not significant" in output.lower() or "ns" in output.lower() or "—" in output


# =============================================================================
# Tests for F3 formatting
# =============================================================================

class TestF3Formatting:
    """Tests for F3 result formatting."""

    def test_format_f3_result_contains_populations(self):
        """F3 output should show target and source populations."""
        from popmix.cli.rich_output import format_f3_result

        result = MockF3Result(target="French", source1="Yamnaya", source2="Anatolia_N")
        output = format_f3_result(result)

        assert "French" in output
        assert "Yamnaya" in output
        assert "Anatolia_N" in output

    def test_format_f3_negative_indicates_admixture(self):
        """Negative F3 with significant Z should indicate admixture."""
        from popmix.cli.rich_output import format_f3_result

        result = MockF3Result(estimate=-0.001, z_score=-10.0)
        output = format_f3_result(result)

        assert "admix" in output.lower()


# =============================================================================
# Tests for F4 formatting
# =============================================================================

class TestF4Formatting:
    """Tests for F4 result formatting."""

    def test_format_f4_result_contains_all_populations(self):
        """F4 output should show all four populations."""
        from popmix.cli.rich_output import format_f4_result

        result = MockF4Result(pop1="A", pop2="B", pop3="C", pop4="D")
        output = format_f4_result(result)

        assert "A" in output
        assert "B" in output
        assert "C" in output
        assert "D" in output

    def test_format_f4_significant_indicates_gene_flow(self):
        """Significant F4 should indicate deviation from tree."""
        from popmix.cli.rich_output import format_f4_result

        result = MockF4Result(z_score=5.0)
        output = format_f4_result(result)

        assert "gene flow" in output.lower() or "deviation" in output.lower() or "significant" in output.lower()


# =============================================================================
# Tests for qpAdm formatting
# =============================================================================

class TestQpAdmFormatting:
    """Tests for qpAdm result formatting."""

    def test_format_qpadm_result_contains_target(self):
        """qpAdm output should show target population."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(target="French")
        output = format_qpadm_result(result)

        assert "French" in output

    def test_format_qpadm_result_contains_sources(self):
        """qpAdm output should show all source populations."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(sources=["Yamnaya", "Anatolia_N", "WHG"])
        output = format_qpadm_result(result)

        assert "Yamnaya" in output
        assert "Anatolia_N" in output
        assert "WHG" in output

    def test_format_qpadm_result_contains_weights_as_percentages(self):
        """qpAdm output should show weights as percentages."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(weights=[0.45, 0.55])
        output = format_qpadm_result(result)

        assert "45" in output  # 45%
        assert "55" in output  # 55%

    def test_format_qpadm_result_contains_confidence_intervals(self):
        """qpAdm output should show confidence intervals for weights."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(weights=[0.45, 0.55], std_errs=[0.036, 0.036])
        output = format_qpadm_result(result)

        # Should have CI brackets
        assert "[" in output and "]" in output

    def test_format_qpadm_result_contains_model_fit(self):
        """qpAdm output should show chi-squared, df, and p-value."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(chi_sq=4.23, df=3, p_value=0.237)
        output = format_qpadm_result(result)

        assert "4.23" in output or "χ²" in output or "chi" in output.lower()
        assert "0.237" in output or "0.24" in output

    def test_format_qpadm_good_fit_indicated(self):
        """qpAdm with p > 0.05 should indicate good fit."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(p_value=0.237)
        output = format_qpadm_result(result)

        assert "pass" in output.lower() or "good" in output.lower() or "✓" in output

    def test_format_qpadm_poor_fit_indicated(self):
        """qpAdm with p < 0.05 should indicate poor fit."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(p_value=0.001)
        output = format_qpadm_result(result)

        assert "fail" in output.lower() or "poor" in output.lower() or "✗" in output or "reject" in output.lower()

    def test_format_qpadm_feasibility_shown(self):
        """qpAdm should show feasibility status."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(weights=[0.45, 0.55])
        output = format_qpadm_result(result)

        assert "feasib" in output.lower() or "valid" in output.lower()

    def test_format_qpadm_infeasible_weights_flagged(self):
        """qpAdm with weights outside [0,1] should be flagged."""
        from popmix.cli.rich_output import format_qpadm_result

        result = MockQpAdmResult(weights=[-0.1, 1.1])
        output = format_qpadm_result(result)

        assert "infeasib" in output.lower() or "invalid" in output.lower() or "✗" in output


# =============================================================================
# Tests for qpWave formatting
# =============================================================================

class TestQpWaveFormatting:
    """Tests for qpWave result formatting."""

    def test_format_qpwave_result_contains_populations(self):
        """qpWave output should show left and right populations."""
        from popmix.cli.rich_output import format_qpwave_result

        result = MockQpWaveResult(
            left_pops=["French", "Spanish"],
            right_pops=["Mbuti", "Han"],
        )
        output = format_qpwave_result(result)

        assert "French" in output
        assert "Spanish" in output

    def test_format_qpwave_result_contains_rank_tests(self):
        """qpWave output should show rank test results."""
        from popmix.cli.rich_output import format_qpwave_result

        result = MockQpWaveResult()
        output = format_qpwave_result(result)

        assert "rank" in output.lower() or "Rank" in output

    def test_format_qpwave_result_shows_estimated_rank(self):
        """qpWave output should show estimated rank."""
        from popmix.cli.rich_output import format_qpwave_result

        result = MockQpWaveResult(estimated_rank=1)
        output = format_qpwave_result(result)

        # Should indicate rank or number of streams
        assert "1" in output or "2" in output  # rank 1 = 2 streams


# =============================================================================
# Tests for DATES formatting
# =============================================================================

class TestDatesFormatting:
    """Tests for DATES result formatting."""

    def test_format_dates_result_contains_populations(self):
        """DATES output should show test and source populations."""
        from popmix.cli.rich_output import format_dates_result

        result = MockDatesResult(
            test_pop="French",
            source1="Yamnaya",
            source2="Anatolia_N",
        )
        output = format_dates_result(result)

        assert "French" in output
        assert "Yamnaya" in output
        assert "Anatolia_N" in output

    def test_format_dates_result_contains_generations(self):
        """DATES output should show generations estimate with SE."""
        from popmix.cli.rich_output import format_dates_result

        result = MockDatesResult(generations=150.2, generations_se=12.3)
        output = format_dates_result(result)

        assert "150" in output
        assert "12" in output or "±" in output

    def test_format_dates_result_contains_years_bp(self):
        """DATES output should show years BP estimate."""
        from popmix.cli.rich_output import format_dates_result

        result = MockDatesResult(years_bp=4356.0, years_bp_se=357.0)
        output = format_dates_result(result)

        assert "4356" in output or "4,356" in output

    def test_format_dates_result_contains_fit_quality(self):
        """DATES output should show R² fit quality."""
        from popmix.cli.rich_output import format_dates_result

        result = MockDatesResult(r_squared=0.92)
        output = format_dates_result(result)

        assert "0.92" in output or "R²" in output or "r2" in output.lower()


# =============================================================================
# Tests for significance coloring
# =============================================================================

class TestSignificanceColoring:
    """Tests for Z-score significance coloring."""

    def test_get_significance_color_high_z(self):
        """Z >= 3 should return green for significant."""
        from popmix.cli.rich_output import get_significance_style

        style = get_significance_style(5.0)
        assert "green" in style.lower() or style == "significant"

    def test_get_significance_color_marginal_z(self):
        """2 <= Z < 3 should return yellow for marginal."""
        from popmix.cli.rich_output import get_significance_style

        style = get_significance_style(2.5)
        assert "yellow" in style.lower() or style == "marginal"

    def test_get_significance_color_low_z(self):
        """Z < 2 should return gray for not significant."""
        from popmix.cli.rich_output import get_significance_style

        style = get_significance_style(1.0)
        assert "gray" in style.lower() or "dim" in style.lower() or style == "not_significant"

    def test_get_significance_color_negative_z(self):
        """Negative Z should use absolute value."""
        from popmix.cli.rich_output import get_significance_style

        style_positive = get_significance_style(5.0)
        style_negative = get_significance_style(-5.0)

        # Both should be significant
        assert style_positive == style_negative


# =============================================================================
# Tests for confidence interval calculation
# =============================================================================

class TestConfidenceInterval:
    """Tests for confidence interval calculation."""

    def test_calculate_ci_95_percent(self):
        """95% CI should use 1.96 * SE."""
        from popmix.cli.rich_output import calculate_confidence_interval

        lower, upper = calculate_confidence_interval(0.02, 0.001, confidence=0.95)

        # 0.02 ± 1.96 * 0.001 = [0.01804, 0.02196]
        assert abs(lower - 0.01804) < 0.0001
        assert abs(upper - 0.02196) < 0.0001

    def test_calculate_ci_99_percent(self):
        """99% CI should use 2.576 * SE."""
        from popmix.cli.rich_output import calculate_confidence_interval

        lower, upper = calculate_confidence_interval(0.02, 0.001, confidence=0.99)

        # 0.02 ± 2.576 * 0.001 = [0.017424, 0.022576]
        assert abs(lower - 0.017424) < 0.0001
        assert abs(upper - 0.022576) < 0.0001


# =============================================================================
# Tests for Rich table creation
# =============================================================================

class TestRichTables:
    """Tests for Rich table creation utilities."""

    def test_create_population_table(self):
        """Population table should format population counts."""
        from popmix.cli.rich_output import create_population_table

        populations = [("French.DG", 45), ("Han.DG", 32), ("Yoruba.DG", 28)]
        table = create_population_table(populations)

        # Should return a Rich Table
        from rich.table import Table
        assert isinstance(table, Table)

    def test_create_fstat_table(self):
        """F-stat table should include all columns."""
        from popmix.cli.rich_output import create_fstat_table

        results = [
            {"pop1": "A", "pop2": "B", "estimate": 0.02, "std_err": 0.001, "z_score": 20.0},
        ]
        table = create_fstat_table(results)

        from rich.table import Table
        assert isinstance(table, Table)


# =============================================================================
# Tests for progress bar
# =============================================================================

class TestProgressBar:
    """Tests for progress bar utilities."""

    def test_create_progress_context(self):
        """Progress context should be creatable."""
        from popmix.cli.rich_output import create_progress

        progress = create_progress("Processing SNPs")
        assert progress is not None
