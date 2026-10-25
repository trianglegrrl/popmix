"""Data aggregation for MultiQC popmix module.

Aggregates results from multiple samples for summary tables.
"""

from typing import Any


def aggregate_results(results: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Aggregate results by type.

    Parameters
    ----------
    results : List[Dict[str, Any]]
        List of result dictionaries

    Returns
    -------
    Dict[str, List[Dict[str, Any]]]
        Results grouped by type
    """
    aggregated = {}

    for result in results:
        result_type = result.get("type", "unknown")

        if result_type not in aggregated:
            aggregated[result_type] = []

        aggregated[result_type].append(result)

    return aggregated


def summarize_f2_results(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize F2 results.

    Parameters
    ----------
    results : List[Dict[str, Any]]
        List of F2 results

    Returns
    -------
    Dict[str, Any]
        Summary statistics
    """
    if not results:
        return {}

    estimates = [r.get("estimate", 0) for r in results]
    z_scores = [abs(r.get("z_score", 0)) for r in results]

    return {
        "n_comparisons": len(results),
        "mean_f2": sum(estimates) / len(estimates) if estimates else 0,
        "max_f2": max(estimates) if estimates else 0,
        "min_f2": min(estimates) if estimates else 0,
        "n_significant": sum(1 for z in z_scores if z >= 2),
    }


def summarize_qpadm_results(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize qpAdm results.

    Parameters
    ----------
    results : List[Dict[str, Any]]
        List of qpAdm results

    Returns
    -------
    Dict[str, Any]
        Summary statistics
    """
    if not results:
        return {}

    p_values = [r.get("p_value", 0) for r in results]
    feasible = [r.get("feasible", True) for r in results]

    return {
        "n_models": len(results),
        "n_good_fit": sum(1 for p in p_values if p > 0.05),
        "n_feasible": sum(1 for f in feasible if f),
        "mean_pvalue": sum(p_values) / len(p_values) if p_values else 0,
    }


def summarize_dates_results(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize DATES results.

    Parameters
    ----------
    results : List[Dict[str, Any]]
        List of DATES results

    Returns
    -------
    Dict[str, Any]
        Summary statistics
    """
    if not results:
        return {}

    generations = [r.get("generations", 0) for r in results]
    r_squared = [r.get("r_squared", 0) for r in results]

    return {
        "n_analyses": len(results),
        "mean_generations": sum(generations) / len(generations) if generations else 0,
        "range_generations": (min(generations), max(generations)) if generations else (0, 0),
        "mean_r_squared": sum(r_squared) / len(r_squared) if r_squared else 0,
    }
