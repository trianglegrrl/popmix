"""Rich terminal output formatting for popmix CLI.

This module provides beautiful, informative terminal output using the Rich library.
All formatters include:
- Colored significance indicators
- Confidence intervals
- Structured tables and panels
"""

from typing import Any, List, Optional, Tuple, Union
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich.text import Text
from rich import box
import math

# Create a shared console instance
console = Console()


# =============================================================================
# Confidence Interval Calculation
# =============================================================================

def calculate_confidence_interval(
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
        Lower and upper bounds of the confidence interval
    """
    # Z-scores for common confidence levels
    z_scores = {
        0.90: 1.645,
        0.95: 1.96,
        0.99: 2.576,
    }
    
    z = z_scores.get(confidence, 1.96)
    margin = z * std_err
    
    return (estimate - margin, estimate + margin)


# =============================================================================
# Significance Styling
# =============================================================================

def get_significance_style(z_score: float) -> str:
    """Get Rich style based on Z-score significance.
    
    Parameters
    ----------
    z_score : float
        Z-score value
    
    Returns
    -------
    str
        Rich style string
    """
    abs_z = abs(z_score)
    
    if abs_z >= 3:
        return "bold green"
    elif abs_z >= 2:
        return "yellow"
    else:
        return "dim"


def get_significance_label(z_score: float) -> str:
    """Get significance label based on Z-score.
    
    Parameters
    ----------
    z_score : float
        Z-score value
    
    Returns
    -------
    str
        Significance label
    """
    abs_z = abs(z_score)
    
    if abs_z >= 3:
        return "✓ Significant"
    elif abs_z >= 2:
        return "~ Marginal"
    else:
        return "— Not significant"


def get_pvalue_style(p_value: float) -> str:
    """Get Rich style based on p-value.
    
    Parameters
    ----------
    p_value : float
        P-value
    
    Returns
    -------
    str
        Rich style string
    """
    if p_value > 0.05:
        return "bold green"
    elif p_value > 0.01:
        return "yellow"
    else:
        return "bold red"


# =============================================================================
# F2 Formatting
# =============================================================================

def format_f2_result(result: Any, show_ci: bool = True) -> str:
    """Format F2 result for Rich display.
    
    Parameters
    ----------
    result : MockF2Result or similar
        F2 result object with estimate, std_err, z_score, n_snps
    show_ci : bool
        Whether to show confidence interval
    
    Returns
    -------
    str
        Formatted string for display
    """
    lines = []
    
    # Header
    lines.append(f"f2({result.pop1}, {result.pop2})")
    lines.append("─" * 50)
    
    # Estimate with CI
    if show_ci and hasattr(result, 'std_err'):
        lower, upper = calculate_confidence_interval(result.estimate, result.std_err)
        lines.append(f"  Estimate:  {result.estimate:.6f}  [{lower:.6f}, {upper:.6f}]")
    else:
        lines.append(f"  Estimate:  {result.estimate:.6f}")
    
    # Standard error
    lines.append(f"  Std Error: {result.std_err:.6f}")
    
    # Z-score with significance
    sig_label = get_significance_label(result.z_score)
    lines.append(f"  Z-score:   {result.z_score:.3f}  {sig_label}")
    
    # SNP count
    lines.append(f"  N SNPs:    {int(result.n_snps):,}")
    
    return "\n".join(lines)


def format_f2_result_rich(result: Any, show_ci: bool = True) -> Panel:
    """Format F2 result as a Rich Panel.
    
    Parameters
    ----------
    result : MockF2Result or similar
        F2 result object
    show_ci : bool
        Whether to show confidence interval
    
    Returns
    -------
    Panel
        Rich Panel object
    """
    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Label", style="bold")
    table.add_column("Value")
    table.add_column("Extra", style="dim")
    
    # Estimate with CI
    if show_ci and hasattr(result, 'std_err'):
        lower, upper = calculate_confidence_interval(result.estimate, result.std_err)
        ci_str = f"[{lower:.6f}, {upper:.6f}]"
    else:
        ci_str = ""
    
    table.add_row("Estimate", f"{result.estimate:.6f}", ci_str)
    table.add_row("Std Error", f"{result.std_err:.6f}", "")
    
    # Z-score with styling
    sig_style = get_significance_style(result.z_score)
    sig_label = get_significance_label(result.z_score)
    z_text = Text(f"{result.z_score:.3f}")
    z_text.stylize(sig_style)
    table.add_row("Z-score", z_text, sig_label)
    
    table.add_row("N SNPs", f"{int(result.n_snps):,}", "")
    
    return Panel(
        table,
        title=f"f2({result.pop1}, {result.pop2})",
        border_style="blue",
    )


# =============================================================================
# F3 Formatting
# =============================================================================

def format_f3_result(result: Any, show_ci: bool = True) -> str:
    """Format F3 result for Rich display.
    
    Parameters
    ----------
    result : MockF3Result or similar
        F3 result object
    show_ci : bool
        Whether to show confidence interval
    
    Returns
    -------
    str
        Formatted string for display
    """
    lines = []
    
    # Header
    lines.append(f"f3({result.target}; {result.source1}, {result.source2})")
    lines.append("─" * 50)
    
    # Estimate with CI
    if show_ci and hasattr(result, 'std_err'):
        lower, upper = calculate_confidence_interval(result.estimate, result.std_err)
        lines.append(f"  Estimate:  {result.estimate:.6f}  [{lower:.6f}, {upper:.6f}]")
    else:
        lines.append(f"  Estimate:  {result.estimate:.6f}")
    
    # Standard error
    lines.append(f"  Std Error: {result.std_err:.6f}")
    
    # Z-score with significance
    sig_label = get_significance_label(result.z_score)
    lines.append(f"  Z-score:   {result.z_score:.3f}  {sig_label}")
    
    # Block count
    lines.append(f"  N blocks:  {int(result.n_blocks):,}")
    
    # Admixture interpretation
    if result.estimate < 0 and abs(result.z_score) >= 3:
        lines.append("")
        lines.append("  => Significant admixture signal detected")
    
    return "\n".join(lines)


# =============================================================================
# F4 Formatting
# =============================================================================

def format_f4_result(result: Any, show_ci: bool = True) -> str:
    """Format F4 result for Rich display.
    
    Parameters
    ----------
    result : MockF4Result or similar
        F4 result object
    show_ci : bool
        Whether to show confidence interval
    
    Returns
    -------
    str
        Formatted string for display
    """
    lines = []
    
    # Header
    lines.append(f"f4({result.pop1}, {result.pop2}; {result.pop3}, {result.pop4})")
    lines.append("─" * 50)
    
    # Estimate with CI
    if show_ci and hasattr(result, 'std_err'):
        lower, upper = calculate_confidence_interval(result.estimate, result.std_err)
        lines.append(f"  Estimate:  {result.estimate:.6f}  [{lower:.6f}, {upper:.6f}]")
    else:
        lines.append(f"  Estimate:  {result.estimate:.6f}")
    
    # Standard error
    lines.append(f"  Std Error: {result.std_err:.6f}")
    
    # Z-score with significance
    sig_label = get_significance_label(result.z_score)
    lines.append(f"  Z-score:   {result.z_score:.3f}  {sig_label}")
    
    # Block count
    lines.append(f"  N blocks:  {int(result.n_blocks):,}")
    
    # Gene flow interpretation
    if abs(result.z_score) >= 3:
        lines.append("")
        lines.append("  => Significant deviation from tree-like evolution (gene flow)")
    
    return "\n".join(lines)


# =============================================================================
# qpAdm Formatting
# =============================================================================

def format_qpadm_result(result: Any, show_ci: bool = True) -> str:
    """Format qpAdm result for Rich display.
    
    Parameters
    ----------
    result : MockQpAdmResult or similar
        qpAdm result object
    show_ci : bool
        Whether to show confidence intervals
    
    Returns
    -------
    str
        Formatted string for display
    """
    lines = []
    
    # Header
    sources_str = " + ".join(result.sources)
    lines.append(f"qpAdm: {result.target} = {sources_str}")
    lines.append("═" * 60)
    
    # Weights table header
    if show_ci and hasattr(result, 'std_errs') and result.std_errs:
        lines.append(f"  {'Source':<20} {'Weight':>10} {'95% CI':>20} {'SE':>10}")
        lines.append("  " + "─" * 56)
        
        for source, weight, se in zip(result.sources, result.weights, result.std_errs):
            lower, upper = calculate_confidence_interval(weight, se)
            pct = weight * 100
            lower_pct = lower * 100
            upper_pct = upper * 100
            ci_str = f"[{lower_pct:.1f}%, {upper_pct:.1f}%]"
            lines.append(f"  {source:<20} {pct:>9.1f}% {ci_str:>20} ±{se*100:.1f}%")
    else:
        lines.append(f"  {'Source':<20} {'Weight':>10}")
        lines.append("  " + "─" * 30)
        
        for source, weight in zip(result.sources, result.weights):
            pct = weight * 100
            lines.append(f"  {source:<20} {pct:>9.1f}%")
    
    lines.append("═" * 60)
    
    # Model fit
    p_status = "✓ PASS" if result.p_value > 0.05 else "✗ FAIL (poor fit)"
    lines.append(f"  Model Fit:    χ²={result.chi_sq:.4f}  df={result.df}  p={result.p_value:.4f}  {p_status}")
    
    # Feasibility
    feasible = result.is_feasible() if hasattr(result, 'is_feasible') else all(0 <= w <= 1 for w in result.weights)
    feas_status = "✓ Valid" if feasible else "✗ Invalid (weights outside [0,1])"
    lines.append(f"  Feasibility:  {feas_status}")
    
    return "\n".join(lines)


def format_qpadm_result_rich(result: Any, show_ci: bool = True) -> Panel:
    """Format qpAdm result as a Rich Panel with table.
    
    Parameters
    ----------
    result : MockQpAdmResult or similar
        qpAdm result object
    show_ci : bool
        Whether to show confidence intervals
    
    Returns
    -------
    Panel
        Rich Panel object
    """
    # Create weights table
    table = Table(box=box.SIMPLE)
    table.add_column("Source", style="bold cyan")
    table.add_column("Weight", justify="right")
    if show_ci and hasattr(result, 'std_errs') and result.std_errs:
        table.add_column("95% CI", justify="right")
        table.add_column("SE", justify="right", style="dim")
    
    for i, source in enumerate(result.sources):
        weight = result.weights[i]
        pct = f"{weight * 100:.1f}%"
        
        if show_ci and hasattr(result, 'std_errs') and result.std_errs:
            se = result.std_errs[i]
            lower, upper = calculate_confidence_interval(weight, se)
            ci_str = f"[{lower*100:.1f}%, {upper*100:.1f}%]"
            se_str = f"±{se*100:.1f}%"
            table.add_row(source, pct, ci_str, se_str)
        else:
            table.add_row(source, pct)
    
    # Model fit status
    p_style = get_pvalue_style(result.p_value)
    p_status = "PASS" if result.p_value > 0.05 else "FAIL"
    
    feasible = result.is_feasible() if hasattr(result, 'is_feasible') else all(0 <= w <= 1 for w in result.weights)
    feas_text = "✓ Valid" if feasible else "✗ Invalid"
    
    # Build content
    from rich.console import Group
    
    fit_text = Text()
    fit_text.append("Model: ", style="bold")
    fit_text.append(f"χ²={result.chi_sq:.2f} df={result.df} ", style="dim")
    fit_text.append(f"p={result.p_value:.4f} ", style=p_style)
    fit_text.append(f"[{p_status}]", style=p_style)
    
    feas_line = Text()
    feas_line.append("Feasibility: ", style="bold")
    feas_line.append(feas_text, style="green" if feasible else "red")
    
    content = Group(table, Text(""), fit_text, feas_line)
    
    sources_str = " + ".join(result.sources)
    return Panel(
        content,
        title=f"[bold]qpAdm: {result.target} = {sources_str}[/bold]",
        border_style="green" if result.p_value > 0.05 else "red",
    )


# =============================================================================
# qpWave Formatting
# =============================================================================

def format_qpwave_result(result: Any) -> str:
    """Format qpWave result for Rich display.
    
    Parameters
    ----------
    result : MockQpWaveResult or similar
        qpWave result object
    
    Returns
    -------
    str
        Formatted string for display
    """
    lines = []
    
    # Header
    left_str = ", ".join(result.left_pops)
    right_str = ", ".join(result.right_pops)
    lines.append(f"qpWave: {left_str}")
    lines.append(f"Outgroups: {right_str}")
    lines.append("═" * 50)
    
    # Rank tests table
    lines.append("")
    lines.append(f"  {'Rank':>6} {'Chi-sq':>10} {'df':>6} {'P-value':>12} {'Status':>8}")
    lines.append("  " + "─" * 44)
    
    for test in result.rank_tests:
        sig = "*" if test.p_value < 0.05 else " "
        status = "Reject" if test.p_value < 0.05 else "Accept"
        lines.append(f"  {test.rank:>6} {test.chi_sq:>10.4f} {test.df:>6} {test.p_value:>12.6f} {status:>8}{sig}")
    
    lines.append("")
    lines.append(f"  Estimated rank: {result.estimated_rank}")
    
    # Interpretation
    if hasattr(result, 'is_clade') and result.is_clade():
        lines.append("  => Populations form a clade (rank 0 not rejected)")
    else:
        n_streams = result.estimated_rank + 1
        lines.append(f"  => {n_streams} independent ancestry streams")
    
    return "\n".join(lines)


# =============================================================================
# DATES Formatting
# =============================================================================

def format_dates_result(result: Any) -> str:
    """Format DATES result for Rich display.
    
    Parameters
    ----------
    result : MockDatesResult or similar
        DATES result object
    
    Returns
    -------
    str
        Formatted string for display
    """
    lines = []
    
    # Header
    lines.append(f"DATES: {result.test_pop} = mix of {result.source1} and {result.source2}")
    lines.append("═" * 60)
    
    # Timing estimates
    lines.append("")
    lines.append(f"  Generations:   {result.generations:.1f} ± {result.generations_se:.1f}")
    
    # Format years with thousands separator
    years_str = f"{int(result.years_bp):,}"
    years_se_str = f"{int(result.years_bp_se):,}"
    lines.append(f"  Years BP:      {years_str} ± {years_se_str} (29 years/gen)")
    
    # Fit quality
    lines.append("")
    lines.append(f"  Amplitude:     {result.amplitude:.6f}")
    lines.append(f"  R² fit:        {result.r_squared:.4f}")
    lines.append(f"  N SNPs:        {result.n_snps:,}")
    
    return "\n".join(lines)


# =============================================================================
# Table Utilities
# =============================================================================

def create_population_table(populations: List[Tuple[str, int]]) -> Table:
    """Create a Rich table for population listing.
    
    Parameters
    ----------
    populations : List[Tuple[str, int]]
        List of (population_name, sample_count) tuples
    
    Returns
    -------
    Table
        Rich Table object
    """
    table = Table(box=box.SIMPLE)
    table.add_column("Population", style="cyan")
    table.add_column("Samples", justify="right", style="green")
    
    for pop, count in populations:
        table.add_row(pop, str(count))
    
    return table


def create_fstat_table(results: List[dict]) -> Table:
    """Create a Rich table for F-statistics results.
    
    Parameters
    ----------
    results : List[dict]
        List of result dictionaries with pop1, pop2, estimate, std_err, z_score
    
    Returns
    -------
    Table
        Rich Table object
    """
    table = Table(box=box.SIMPLE)
    table.add_column("Pop1", style="cyan")
    table.add_column("Pop2", style="cyan")
    table.add_column("Estimate", justify="right")
    table.add_column("Std Err", justify="right", style="dim")
    table.add_column("Z-score", justify="right")
    table.add_column("Significance")
    
    for r in results:
        sig_style = get_significance_style(r["z_score"])
        sig_label = get_significance_label(r["z_score"])
        
        z_text = Text(f"{r['z_score']:.2f}")
        z_text.stylize(sig_style)
        
        table.add_row(
            r["pop1"],
            r["pop2"],
            f"{r['estimate']:.6f}",
            f"{r['std_err']:.6f}",
            z_text,
            sig_label,
        )
    
    return table


# =============================================================================
# Progress Bar Utilities
# =============================================================================

def create_progress(description: str = "Processing") -> Progress:
    """Create a Rich progress bar.
    
    Parameters
    ----------
    description : str
        Task description
    
    Returns
    -------
    Progress
        Rich Progress object
    """
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    )


def print_panel(content: str, title: str = "", style: str = "blue") -> None:
    """Print content in a Rich panel.
    
    Parameters
    ----------
    content : str
        Panel content
    title : str
        Panel title
    style : str
        Border style
    """
    console.print(Panel(content, title=title, border_style=style))
