"""General Statistics columns for MultiQC popmix module.

Defines the columns that appear in the MultiQC General Statistics table.
"""

from typing import Any, Dict


def get_general_stats_headers() -> Dict[str, Dict[str, Any]]:
    """Get General Statistics column headers.
    
    Returns
    -------
    Dict[str, Dict[str, Any]]
        Column header definitions
    """
    return {
        "popmix_pvalue": {
            "title": "P-value",
            "description": "Model fit p-value (qpAdm/qpWave)",
            "min": 0,
            "max": 1,
            "scale": "RdYlGn",
            "format": "{:,.3f}",
        },
        "popmix_generations": {
            "title": "Generations",
            "description": "Admixture date (DATES)",
            "min": 0,
            "scale": "Blues",
            "format": "{:,.0f}",
        },
        "popmix_years_bp": {
            "title": "Years BP",
            "description": "Years before present (DATES)",
            "min": 0,
            "scale": "Purples",
            "format": "{:,.0f}",
        },
        "popmix_r_squared": {
            "title": "R²",
            "description": "Model fit R-squared",
            "min": 0,
            "max": 1,
            "scale": "Greens",
            "format": "{:,.2f}",
        },
        "popmix_n_sources": {
            "title": "# Sources",
            "description": "Number of source populations",
            "min": 0,
            "scale": "Greys",
            "format": "{:,.0f}",
        },
    }


def format_general_stats_data(results: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Format results for General Statistics table.
    
    Parameters
    ----------
    results : Dict[str, Dict[str, Any]]
        Results by sample name
    
    Returns
    -------
    Dict[str, Dict[str, Any]]
        Formatted data for MultiQC
    """
    formatted = {}
    
    for sample_name, data in results.items():
        sample_data = {}
        
        # Extract relevant fields
        if "p_value" in data:
            sample_data["popmix_pvalue"] = data["p_value"]
        
        if "generations" in data:
            sample_data["popmix_generations"] = data["generations"]
        
        if "years_bp" in data:
            sample_data["popmix_years_bp"] = data["years_bp"]
        
        if "r_squared" in data:
            sample_data["popmix_r_squared"] = data["r_squared"]
        
        if "sources" in data:
            sample_data["popmix_n_sources"] = len(data["sources"])
        
        if sample_data:
            formatted[sample_name] = sample_data
    
    return formatted
