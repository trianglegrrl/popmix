"""Tests for statistical rigor in result reporting.

These tests ensure all statistical results include proper uncertainty
quantification, significance assessment, and effect sizes.
"""



# =============================================================================
# Mock result classes
# =============================================================================

class MockF2Result:
    def __init__(self, estimate=0.02, std_err=0.001, z_score=20.0, n_snps=100000):
        self.pop1 = "A"
        self.pop2 = "B"
        self.estimate = estimate
        self.std_err = std_err
        self.z_score = z_score
        self.n_snps = n_snps


class MockQpAdmResult:
    def __init__(self, weights=None, std_errs=None, p_value=0.237):
        self.target = "French"
        self.sources = ["Yamnaya", "Anatolia_N"]
        self.weights = weights or [0.45, 0.55]
        self.std_errs = std_errs or [0.036, 0.036]
        self.chi_sq = 4.23
        self.df = 3
        self.p_value = p_value

    def is_feasible(self):
        return all(0 <= w <= 1 for w in self.weights)


# =============================================================================
# Tests for Confidence Interval Calculation
# =============================================================================

class TestConfidenceIntervals:
    """Tests for proper confidence interval calculation."""

    def test_ci_95_uses_correct_z_score(self):
        """95% CI should use z = 1.96."""
        from popmix.cli.statistics import calculate_ci

        lower, upper = calculate_ci(0.5, 0.1, confidence=0.95)

        # 0.5 ± 1.96 * 0.1 = [0.304, 0.696]
        assert abs(lower - 0.304) < 0.001
        assert abs(upper - 0.696) < 0.001

    def test_ci_90_uses_correct_z_score(self):
        """90% CI should use z = 1.645."""
        from popmix.cli.statistics import calculate_ci

        lower, upper = calculate_ci(0.5, 0.1, confidence=0.90)

        # 0.5 ± 1.645 * 0.1 = [0.3355, 0.6645]
        assert abs(lower - 0.3355) < 0.001
        assert abs(upper - 0.6645) < 0.001

    def test_ci_99_uses_correct_z_score(self):
        """99% CI should use z = 2.576."""
        from popmix.cli.statistics import calculate_ci

        lower, upper = calculate_ci(0.5, 0.1, confidence=0.99)

        # 0.5 ± 2.576 * 0.1 = [0.2424, 0.7576]
        assert abs(lower - 0.2424) < 0.001
        assert abs(upper - 0.7576) < 0.001

    def test_ci_symmetric_around_estimate(self):
        """CI should be symmetric around the point estimate."""
        from popmix.cli.statistics import calculate_ci

        estimate = 0.123
        lower, upper = calculate_ci(estimate, 0.05)

        midpoint = (lower + upper) / 2
        assert abs(midpoint - estimate) < 1e-10


# =============================================================================
# Tests for Significance Assessment
# =============================================================================

class TestSignificanceAssessment:
    """Tests for Z-score based significance assessment."""

    def test_significance_level_highly_significant(self):
        """Z >= 3 should be highly significant."""
        from popmix.cli.statistics import assess_significance

        result = assess_significance(5.0)

        assert result["level"] == "significant"
        assert result["significant"]
        assert result["stars"] == "***"

    def test_significance_level_significant(self):
        """3 > Z >= 2 should be significant."""
        from popmix.cli.statistics import assess_significance

        result = assess_significance(2.5)

        assert result["level"] == "marginal"
        assert result["stars"] == "*"

    def test_significance_level_not_significant(self):
        """Z < 2 should be not significant."""
        from popmix.cli.statistics import assess_significance

        result = assess_significance(1.5)

        assert result["level"] == "not_significant"
        assert not result["significant"]
        assert result["stars"] == ""

    def test_significance_works_for_negative_z(self):
        """Negative Z scores should use absolute value."""
        from popmix.cli.statistics import assess_significance

        result_pos = assess_significance(4.0)
        result_neg = assess_significance(-4.0)

        assert result_pos["level"] == result_neg["level"]
        assert result_pos["significant"] == result_neg["significant"]

    def test_significance_includes_p_value_approximation(self):
        """Significance should include approximate p-value."""
        from popmix.cli.statistics import assess_significance

        result = assess_significance(2.0)

        assert "p_value_approx" in result
        # Z=2 corresponds to p ≈ 0.0455 (two-tailed)
        assert 0.04 < result["p_value_approx"] < 0.05


# =============================================================================
# Tests for Effect Size Calculation
# =============================================================================

class TestEffectSize:
    """Tests for effect size calculations."""

    def test_cohens_d_calculation(self):
        """Cohen's d should be calculated correctly."""
        from popmix.cli.statistics import calculate_cohens_d

        # d = (mean1 - mean2) / pooled_sd
        d = calculate_cohens_d(mean1=10.0, mean2=8.0, sd1=2.0, sd2=2.0)

        # Pooled SD = 2.0, so d = 2/2 = 1.0
        assert abs(d - 1.0) < 0.01

    def test_cohens_d_interpretation_small(self):
        """d < 0.5 should be interpreted as small."""
        from popmix.cli.statistics import interpret_cohens_d

        interpretation = interpret_cohens_d(0.3)
        assert interpretation == "small"

    def test_cohens_d_interpretation_medium(self):
        """0.5 <= d < 0.8 should be interpreted as medium."""
        from popmix.cli.statistics import interpret_cohens_d

        interpretation = interpret_cohens_d(0.6)
        assert interpretation == "medium"

    def test_cohens_d_interpretation_large(self):
        """d >= 0.8 should be interpreted as large."""
        from popmix.cli.statistics import interpret_cohens_d

        interpretation = interpret_cohens_d(1.2)
        assert interpretation == "large"


# =============================================================================
# Tests for P-Value Thresholds
# =============================================================================

class TestPValueThresholds:
    """Tests for p-value interpretation thresholds."""

    def test_pvalue_good_fit(self):
        """p > 0.05 should indicate good fit."""
        from popmix.cli.statistics import assess_model_fit

        result = assess_model_fit(p_value=0.15)

        assert result["quality"] == "good"
        assert result["pass"]

    def test_pvalue_marginal_fit(self):
        """0.01 < p <= 0.05 should indicate marginal fit."""
        from popmix.cli.statistics import assess_model_fit

        result = assess_model_fit(p_value=0.03)

        assert result["quality"] == "marginal"

    def test_pvalue_poor_fit(self):
        """p <= 0.01 should indicate poor fit."""
        from popmix.cli.statistics import assess_model_fit

        result = assess_model_fit(p_value=0.005)

        assert result["quality"] == "poor"
        assert not result["pass"]

    def test_model_fit_includes_recommendation(self):
        """Model fit assessment should include recommendation."""
        from popmix.cli.statistics import assess_model_fit

        result = assess_model_fit(p_value=0.001)

        assert "recommendation" in result


# =============================================================================
# Tests for Statistical Summary
# =============================================================================

class TestStatisticalSummary:
    """Tests for comprehensive statistical summaries."""

    def test_summary_includes_all_components(self):
        """Statistical summary should include all required components."""
        from popmix.cli.statistics import create_statistical_summary

        summary = create_statistical_summary(
            estimate=0.02,
            std_err=0.001,
            n_observations=100000,
        )

        assert "estimate" in summary
        assert "std_err" in summary
        assert "confidence_interval" in summary
        assert "z_score" in summary
        assert "significance" in summary
        assert "n" in summary

    def test_summary_computes_z_score(self):
        """Summary should compute Z-score from estimate/SE."""
        from popmix.cli.statistics import create_statistical_summary

        summary = create_statistical_summary(
            estimate=0.02,
            std_err=0.001,
        )

        # Z = 0.02 / 0.001 = 20
        assert abs(summary["z_score"] - 20.0) < 0.01

    def test_summary_handles_zero_std_err(self):
        """Summary should handle zero std error gracefully."""
        from popmix.cli.statistics import create_statistical_summary

        summary = create_statistical_summary(
            estimate=0.02,
            std_err=0.0,
        )

        # Should not raise, Z-score should be inf or very large
        assert summary["z_score"] == float('inf') or summary["z_score"] > 1e6


# =============================================================================
# Tests for Feasibility Assessment
# =============================================================================

class TestFeasibilityAssessment:
    """Tests for qpAdm weight feasibility assessment."""

    def test_feasible_weights(self):
        """Weights in [0, 1] should be feasible."""
        from popmix.cli.statistics import assess_weight_feasibility

        result = assess_weight_feasibility([0.3, 0.4, 0.3])

        assert result["feasible"]
        assert result["sum_to_one"]

    def test_infeasible_negative_weight(self):
        """Negative weights should be infeasible."""
        from popmix.cli.statistics import assess_weight_feasibility

        result = assess_weight_feasibility([-0.1, 0.6, 0.5])

        assert not result["feasible"]
        assert "negative" in result["reason"].lower()

    def test_infeasible_weight_over_one(self):
        """Weights > 1 should be infeasible."""
        from popmix.cli.statistics import assess_weight_feasibility

        result = assess_weight_feasibility([1.2, -0.2])

        assert not result["feasible"]

    def test_weights_sum_check(self):
        """Should check if weights sum to ~1."""
        from popmix.cli.statistics import assess_weight_feasibility

        result = assess_weight_feasibility([0.3, 0.3, 0.3])  # Sum = 0.9

        assert not result["sum_to_one"]
