"""Statistical rigor utilities for popmix results.

This module provides functions for:
- Confidence interval calculation
- Significance assessment with proper thresholds
- Effect size calculation and interpretation
- Model fit assessment
- Comprehensive statistical summaries
"""

import math
from typing import Any, Dict, List, Optional, Tuple


# =============================================================================
# Confidence Interval Calculation
# =============================================================================

# Z-scores for common confidence levels
Z_SCORES = {
    0.90: 1.645,
    0.95: 1.96,
    0.99: 2.576,
}


def calculate_ci(
    estimate: float,
    std_err: float,
    confidence: float = 0.95,
) -> Tuple[float, float]:
    """Calculate confidence interval for an estimate.
    
    Parameters
    ----------
    estimate : float
        Point estimate
    std_err : float
        Standard error of the estimate
    confidence : float
        Confidence level (default 0.95 for 95% CI)
    
    Returns
    -------
    Tuple[float, float]
        (lower, upper) bounds of the confidence interval
    
    Examples
    --------
    >>> lower, upper = calculate_ci(0.5, 0.1, 0.95)
    >>> print(f"95% CI: [{lower:.3f}, {upper:.3f}]")
    95% CI: [0.304, 0.696]
    """
    z = Z_SCORES.get(confidence, 1.96)
    margin = z * std_err
    return (estimate - margin, estimate + margin)


# =============================================================================
# Significance Assessment
# =============================================================================

def assess_significance(z_score: float) -> Dict[str, Any]:
    """Assess statistical significance based on Z-score.
    
    Uses the following thresholds:
    - |Z| >= 3: Highly significant (***)
    - 2 <= |Z| < 3: Marginally significant (*)
    - |Z| < 2: Not significant
    
    Parameters
    ----------
    z_score : float
        Z-score value
    
    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - level: "significant", "marginal", or "not_significant"
        - significant: bool
        - stars: significance stars ("***", "*", or "")
        - description: human-readable description
        - p_value_approx: approximate two-tailed p-value
    """
    abs_z = abs(z_score)
    
    # Approximate p-value using normal distribution
    # Using the approximation: p ≈ 2 * (1 - Φ(|z|))
    # For standard normal, we use a simple approximation
    if abs_z > 8:
        p_approx = 0.0
    else:
        # Approximation of 1 - Φ(z) for standard normal
        # Using Abramowitz and Stegun approximation
        t = 1 / (1 + 0.2316419 * abs_z)
        d = 0.3989423 * math.exp(-abs_z * abs_z / 2)
        p_one_tail = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
        p_approx = 2 * p_one_tail  # Two-tailed
    
    if abs_z >= 3:
        return {
            "level": "significant",
            "significant": True,
            "stars": "***",
            "description": "Highly significant (p < 0.003)",
            "p_value_approx": p_approx,
            "color": "green",
        }
    elif abs_z >= 2:
        return {
            "level": "marginal",
            "significant": True,
            "stars": "*",
            "description": "Marginally significant (p < 0.05)",
            "p_value_approx": p_approx,
            "color": "yellow",
        }
    else:
        return {
            "level": "not_significant",
            "significant": False,
            "stars": "",
            "description": "Not significant (p >= 0.05)",
            "p_value_approx": p_approx,
            "color": "dim",
        }


# =============================================================================
# Effect Size Calculation
# =============================================================================

def calculate_cohens_d(
    mean1: float,
    mean2: float,
    sd1: float,
    sd2: float,
    n1: Optional[int] = None,
    n2: Optional[int] = None,
) -> float:
    """Calculate Cohen's d effect size.
    
    Parameters
    ----------
    mean1 : float
        Mean of group 1
    mean2 : float
        Mean of group 2
    sd1 : float
        Standard deviation of group 1
    sd2 : float
        Standard deviation of group 2
    n1, n2 : int, optional
        Sample sizes (for weighted pooled SD)
    
    Returns
    -------
    float
        Cohen's d effect size
    """
    if n1 is not None and n2 is not None:
        # Weighted pooled standard deviation
        pooled_var = ((n1 - 1) * sd1**2 + (n2 - 1) * sd2**2) / (n1 + n2 - 2)
        pooled_sd = math.sqrt(pooled_var)
    else:
        # Simple pooled SD
        pooled_sd = math.sqrt((sd1**2 + sd2**2) / 2)
    
    if pooled_sd == 0:
        return float('inf') if mean1 != mean2 else 0.0
    
    return (mean1 - mean2) / pooled_sd


def interpret_cohens_d(d: float) -> str:
    """Interpret Cohen's d effect size.
    
    Uses Cohen's (1988) conventions:
    - |d| < 0.5: small
    - 0.5 <= |d| < 0.8: medium
    - |d| >= 0.8: large
    
    Parameters
    ----------
    d : float
        Cohen's d value
    
    Returns
    -------
    str
        Interpretation: "small", "medium", or "large"
    """
    abs_d = abs(d)
    
    if abs_d < 0.5:
        return "small"
    elif abs_d < 0.8:
        return "medium"
    else:
        return "large"


# =============================================================================
# Model Fit Assessment
# =============================================================================

def assess_model_fit(
    p_value: float,
    chi_sq: Optional[float] = None,
    df: Optional[int] = None,
) -> Dict[str, Any]:
    """Assess model fit quality based on p-value.
    
    Uses the following thresholds:
    - p > 0.05: Good fit
    - 0.01 < p <= 0.05: Marginal fit
    - p <= 0.01: Poor fit
    
    Parameters
    ----------
    p_value : float
        P-value from chi-squared test
    chi_sq : float, optional
        Chi-squared statistic
    df : int, optional
        Degrees of freedom
    
    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - quality: "good", "marginal", or "poor"
        - pass: bool
        - recommendation: interpretation and recommendation
    """
    if p_value > 0.05:
        return {
            "quality": "good",
            "pass": True,
            "color": "green",
            "symbol": "✓",
            "recommendation": "Model fits well. Results are reliable.",
            "p_value": p_value,
            "chi_sq": chi_sq,
            "df": df,
        }
    elif p_value > 0.01:
        return {
            "quality": "marginal",
            "pass": True,
            "color": "yellow",
            "symbol": "~",
            "recommendation": "Model fit is marginal. Interpret with caution.",
            "p_value": p_value,
            "chi_sq": chi_sq,
            "df": df,
        }
    else:
        return {
            "quality": "poor",
            "pass": False,
            "color": "red",
            "symbol": "✗",
            "recommendation": "Model fit is poor. Consider alternative models or additional sources.",
            "p_value": p_value,
            "chi_sq": chi_sq,
            "df": df,
        }


# =============================================================================
# Statistical Summary
# =============================================================================

def create_statistical_summary(
    estimate: float,
    std_err: float,
    n_observations: Optional[int] = None,
    confidence: float = 0.95,
) -> Dict[str, Any]:
    """Create comprehensive statistical summary.
    
    Parameters
    ----------
    estimate : float
        Point estimate
    std_err : float
        Standard error
    n_observations : int, optional
        Number of observations/SNPs/blocks
    confidence : float
        Confidence level for CI
    
    Returns
    -------
    Dict[str, Any]
        Comprehensive summary including CI, Z-score, significance
    """
    # Calculate Z-score
    if std_err == 0:
        z_score = float('inf') if estimate != 0 else 0.0
    else:
        z_score = estimate / std_err
    
    # Calculate confidence interval
    ci_lower, ci_upper = calculate_ci(estimate, std_err, confidence)
    
    # Assess significance
    significance = assess_significance(z_score)
    
    return {
        "estimate": estimate,
        "std_err": std_err,
        "z_score": z_score,
        "confidence_interval": (ci_lower, ci_upper),
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "confidence_level": confidence,
        "significance": significance,
        "n": n_observations,
    }


# =============================================================================
# Weight Feasibility Assessment
# =============================================================================

def assess_weight_feasibility(
    weights: List[float],
    tolerance: float = 0.01,
) -> Dict[str, Any]:
    """Assess feasibility of qpAdm weights.
    
    Parameters
    ----------
    weights : List[float]
        Admixture weights
    tolerance : float
        Tolerance for sum-to-one check
    
    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - feasible: bool (all weights in [0, 1])
        - sum_to_one: bool
        - reason: explanation if infeasible
    """
    # Check if all weights are in [0, 1]
    all_in_range = all(0 <= w <= 1 for w in weights)
    
    # Check if weights sum to approximately 1
    weight_sum = sum(weights)
    sums_to_one = abs(weight_sum - 1.0) < tolerance
    
    # Check for specific issues
    has_negative = any(w < 0 for w in weights)
    has_over_one = any(w > 1 for w in weights)
    
    # Determine reason if infeasible
    reason = ""
    if has_negative:
        reason = "Negative weight(s) detected"
    elif has_over_one:
        reason = "Weight(s) exceed 1.0"
    elif not sums_to_one:
        reason = f"Weights sum to {weight_sum:.3f}, not 1.0"
    
    return {
        "feasible": all_in_range,
        "sum_to_one": sums_to_one,
        "all_in_range": all_in_range,
        "weight_sum": weight_sum,
        "reason": reason,
        "has_negative": has_negative,
        "has_over_one": has_over_one,
    }


# =============================================================================
# Interpretation Helpers
# =============================================================================

def interpret_f3_result(estimate: float, z_score: float) -> Dict[str, Any]:
    """Interpret f3 admixture test result.
    
    Parameters
    ----------
    estimate : float
        f3 estimate
    z_score : float
        Z-score
    
    Returns
    -------
    Dict[str, Any]
        Interpretation dictionary
    """
    significance = assess_significance(z_score)
    
    if estimate < 0 and significance["significant"]:
        return {
            "interpretation": "admixture",
            "description": "Significant negative f3 indicates admixture",
            "confidence": "high" if abs(z_score) >= 3 else "moderate",
        }
    elif estimate < 0:
        return {
            "interpretation": "possible_admixture",
            "description": "Negative f3 but not significant",
            "confidence": "low",
        }
    else:
        return {
            "interpretation": "no_admixture",
            "description": "Positive f3 - no evidence of admixture",
            "confidence": "high" if significance["significant"] else "moderate",
        }


def interpret_f4_result(estimate: float, z_score: float) -> Dict[str, Any]:
    """Interpret f4 / D-statistic result.
    
    Parameters
    ----------
    estimate : float
        f4 estimate
    z_score : float
        Z-score
    
    Returns
    -------
    Dict[str, Any]
        Interpretation dictionary
    """
    significance = assess_significance(z_score)
    
    if significance["significant"]:
        direction = "positive" if estimate > 0 else "negative"
        return {
            "interpretation": "gene_flow",
            "description": f"Significant {direction} f4 indicates gene flow / deviation from tree",
            "direction": direction,
            "confidence": "high" if abs(z_score) >= 3 else "moderate",
        }
    else:
        return {
            "interpretation": "tree_like",
            "description": "No significant deviation from tree-like evolution",
            "confidence": "moderate",
        }
