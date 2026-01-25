"""Batch analysis mode for popmix.

This module provides batch processing for multi-target analyses:
- qpAdm batch: Run qpAdm for multiple targets
- DATES batch: Run DATES for multiple targets

Includes progress tracking, result collection, and summary reports.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Import YAML conditionally
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


# =============================================================================
# Configuration Parsing
# =============================================================================

def parse_populations_file(filepath: str) -> List[str]:
    """Parse population names from a text file.
    
    Parameters
    ----------
    filepath : str
        Path to text file with one population per line
    
    Returns
    -------
    List[str]
        List of population names
    """
    populations = []
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue
            
            populations.append(line)
    
    return populations


def parse_batch_config(filepath: str) -> Dict[str, Any]:
    """Parse batch configuration from YAML file.
    
    Parameters
    ----------
    filepath : str
        Path to YAML config file
    
    Returns
    -------
    Dict[str, Any]
        Configuration dictionary
    """
    if not YAML_AVAILABLE:
        raise ImportError("PyYAML required for config files. Install with: pip install pyyaml")
    
    with open(filepath) as f:
        config = yaml.safe_load(f)
    
    return config


# =============================================================================
# Job Generation
# =============================================================================

def generate_qpadm_jobs(
    targets: List[str],
    sources: List[str],
    outgroups: List[str],
) -> List[Dict[str, Any]]:
    """Generate qpAdm analysis jobs.
    
    Parameters
    ----------
    targets : List[str]
        Target populations to model
    sources : List[str]
        Source populations
    outgroups : List[str]
        Outgroup populations
    
    Returns
    -------
    List[Dict[str, Any]]
        List of job specifications
    """
    jobs = []
    
    for target in targets:
        jobs.append({
            "type": "qpadm",
            "target": target,
            "sources": sources,
            "outgroups": outgroups,
        })
    
    return jobs


def generate_dates_jobs(
    targets: List[str],
    source1: str,
    source2: str,
) -> List[Dict[str, Any]]:
    """Generate DATES analysis jobs.
    
    Parameters
    ----------
    targets : List[str]
        Test (admixed) populations
    source1 : str
        First source population
    source2 : str
        Second source population
    
    Returns
    -------
    List[Dict[str, Any]]
        List of job specifications
    """
    jobs = []
    
    for target in targets:
        jobs.append({
            "type": "dates",
            "test_pop": target,
            "source1": source1,
            "source2": source2,
        })
    
    return jobs


# =============================================================================
# Batch Runner
# =============================================================================

class BatchRunner:
    """Tracks progress of batch analysis jobs.
    
    Attributes
    ----------
    total : int
        Total number of jobs
    completed : int
        Number of completed jobs
    """
    
    def __init__(self, n_jobs: int):
        """Initialize batch runner.
        
        Parameters
        ----------
        n_jobs : int
            Total number of jobs to run
        """
        self.total = n_jobs
        self.completed = 0
        self.failed = 0
    
    def mark_completed(self) -> None:
        """Mark a job as completed."""
        self.completed += 1
    
    def mark_failed(self) -> None:
        """Mark a job as failed."""
        self.failed += 1
        self.completed += 1
    
    @property
    def progress_percent(self) -> float:
        """Get progress as percentage."""
        if self.total == 0:
            return 100.0
        return (self.completed / self.total) * 100
    
    @property
    def remaining(self) -> int:
        """Get number of remaining jobs."""
        return self.total - self.completed


# =============================================================================
# Result Collection
# =============================================================================

class ResultCollector:
    """Collects results from batch analyses.
    
    Attributes
    ----------
    results : List[Dict[str, Any]]
        Collected results
    failures : List[Dict[str, Any]]
        Failed analyses
    """
    
    def __init__(self):
        """Initialize result collector."""
        self._results: List[Dict[str, Any]] = []
        self._failures: List[Dict[str, Any]] = []
    
    def add_result(self, target: str, result: Dict[str, Any]) -> None:
        """Add a successful result.
        
        Parameters
        ----------
        target : str
            Target population name
        result : Dict[str, Any]
            Analysis result
        """
        self._results.append({
            "target": target,
            **result,
        })
    
    def add_failure(self, target: str, error: str) -> None:
        """Add a failed analysis.
        
        Parameters
        ----------
        target : str
            Target population name
        error : str
            Error message
        """
        self._failures.append({
            "target": target,
            "error": error,
        })
    
    def get_results(self) -> List[Dict[str, Any]]:
        """Get all successful results."""
        return self._results
    
    def get_failures(self) -> List[Dict[str, Any]]:
        """Get all failures."""
        return self._failures
    
    @property
    def n_success(self) -> int:
        """Number of successful results."""
        return len(self._results)
    
    @property
    def n_failures(self) -> int:
        """Number of failures."""
        return len(self._failures)


# =============================================================================
# Output Generation
# =============================================================================

def write_batch_results(
    results: List[Dict[str, Any]],
    output_path: str,
    format: str = "json",
) -> None:
    """Write all batch results to a single file.
    
    Parameters
    ----------
    results : List[Dict[str, Any]]
        List of results
    output_path : str
        Output file path
    format : str
        Output format (json, yaml)
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    if format == "json":
        with open(path, "w") as f:
            json.dump(results, f, indent=2)
    elif format == "yaml":
        if not YAML_AVAILABLE:
            raise ImportError("PyYAML required for YAML output")
        with open(path, "w") as f:
            yaml.dump(results, f, default_flow_style=False)
    else:
        raise ValueError(f"Unknown format: {format}")


def write_batch_results_individual(
    results: List[Dict[str, Any]],
    output_dir: str,
    format: str = "json",
    suffix: str = "qpadm",
) -> None:
    """Write individual result files per target.
    
    Parameters
    ----------
    results : List[Dict[str, Any]]
        List of results
    output_dir : str
        Output directory
    format : str
        Output format
    suffix : str
        File suffix (qpadm, dates, etc.)
    """
    dir_path = Path(output_dir)
    dir_path.mkdir(parents=True, exist_ok=True)
    
    for result in results:
        target = result.get("target", "unknown")
        
        if format == "json":
            filename = f"{target}_{suffix}.json"
            filepath = dir_path / filename
            with open(filepath, "w") as f:
                json.dump(result, f, indent=2)
        elif format == "yaml":
            if not YAML_AVAILABLE:
                raise ImportError("PyYAML required for YAML output")
            filename = f"{target}_{suffix}.yaml"
            filepath = dir_path / filename
            with open(filepath, "w") as f:
                yaml.dump(result, f, default_flow_style=False)


# =============================================================================
# Summary Report Generation
# =============================================================================

def generate_batch_summary(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate summary statistics for batch results.
    
    Parameters
    ----------
    results : List[Dict[str, Any]]
        List of analysis results
    
    Returns
    -------
    Dict[str, Any]
        Summary statistics
    """
    if not results:
        return {
            "n_total": 0,
            "n_good_fit": 0,
            "n_feasible": 0,
            "best_models": [],
        }
    
    # Count statistics
    n_total = len(results)
    n_good_fit = sum(1 for r in results if r.get("p_value", 0) > 0.05)
    n_feasible = sum(1 for r in results if r.get("feasible", True))
    
    # Sort by p-value (best fits first)
    sorted_results = sorted(
        results,
        key=lambda r: r.get("p_value", 0),
        reverse=True,
    )
    
    # Get top 5 best models
    best_models = sorted_results[:5]
    
    # Calculate summary statistics
    p_values = [r.get("p_value", 0) for r in results if "p_value" in r]
    
    summary = {
        "n_total": n_total,
        "n_good_fit": n_good_fit,
        "n_feasible": n_feasible,
        "best_models": best_models,
        "mean_pvalue": sum(p_values) / len(p_values) if p_values else 0,
        "max_pvalue": max(p_values) if p_values else 0,
        "min_pvalue": min(p_values) if p_values else 0,
    }
    
    return summary


def format_batch_summary(summary: Dict[str, Any]) -> str:
    """Format batch summary for display.
    
    Parameters
    ----------
    summary : Dict[str, Any]
        Summary statistics
    
    Returns
    -------
    str
        Formatted summary string
    """
    lines = [
        "=" * 60,
        "BATCH ANALYSIS SUMMARY",
        "=" * 60,
        "",
        f"Total analyses: {summary['n_total']}",
        f"Good fit (p > 0.05): {summary['n_good_fit']}",
        f"Feasible models: {summary['n_feasible']}",
        "",
    ]
    
    if summary.get("best_models"):
        lines.append("Best fitting models:")
        for i, model in enumerate(summary["best_models"][:5], 1):
            target = model.get("target", "Unknown")
            p_value = model.get("p_value", 0)
            lines.append(f"  {i}. {target}: p = {p_value:.4f}")
    
    lines.append("")
    lines.append("=" * 60)
    
    return "\n".join(lines)
