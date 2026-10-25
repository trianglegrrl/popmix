"""File parsers for MultiQC popmix module.

Parses popmix JSON output files into structured data.
"""

import json
from pathlib import Path
from typing import Any


def parse_popmix_json(filepath: str) -> dict[str, Any]:
    """Parse a popmix JSON output file.

    Parameters
    ----------
    filepath : str
        Path to JSON file

    Returns
    -------
    Dict[str, Any]
        Parsed data
    """
    with open(filepath) as f:
        data = json.load(f)

    # Infer type from filename if not in data
    if "type" not in data:
        filename = Path(filepath).stem.lower()

        if "qpadm" in filename:
            data["type"] = "qpadm"
        elif "qpwave" in filename:
            data["type"] = "qpwave"
        elif "dates" in filename:
            data["type"] = "dates"
        elif "fst" in filename:
            data["type"] = "fst"
        elif "f4" in filename:
            data["type"] = "f4"
        elif "f3" in filename:
            data["type"] = "f3"
        elif "f2" in filename:
            data["type"] = "f2"
        else:
            # Try to infer from content
            if "target" in data and "sources" in data:
                data["type"] = "qpadm"
            elif "rank_tests" in data:
                data["type"] = "qpwave"
            elif "generations" in data:
                data["type"] = "dates"
            elif "pop1" in data and "pop2" in data:
                data["type"] = "f2"
            else:
                data["type"] = "unknown"

    return data


def parse_f2_result(data: dict[str, Any]) -> dict[str, Any]:
    """Parse F2 result data.

    Parameters
    ----------
    data : Dict[str, Any]
        Raw JSON data

    Returns
    -------
    Dict[str, Any]
        Structured F2 result
    """
    return {
        "type": "f2",
        "pop1": data.get("pop1", ""),
        "pop2": data.get("pop2", ""),
        "estimate": data.get("estimate", 0),
        "std_err": data.get("std_err", 0),
        "z_score": data.get("z_score", 0),
        "n_snps": data.get("n_snps", 0),
    }


def parse_qpadm_result(data: dict[str, Any]) -> dict[str, Any]:
    """Parse qpAdm result data.

    Parameters
    ----------
    data : Dict[str, Any]
        Raw JSON data

    Returns
    -------
    Dict[str, Any]
        Structured qpAdm result
    """
    return {
        "type": "qpadm",
        "target": data.get("target", ""),
        "sources": data.get("sources", []),
        "weights": data.get("weights", {}),
        "chi_sq": data.get("chi_sq", 0),
        "df": data.get("df", 0),
        "p_value": data.get("p_value", 0),
        "feasible": data.get("feasible", True),
    }


def parse_dates_result(data: dict[str, Any]) -> dict[str, Any]:
    """Parse DATES result data.

    Parameters
    ----------
    data : Dict[str, Any]
        Raw JSON data

    Returns
    -------
    Dict[str, Any]
        Structured DATES result
    """
    return {
        "type": "dates",
        "test_pop": data.get("test_pop", ""),
        "source1": data.get("source1", ""),
        "source2": data.get("source2", ""),
        "generations": data.get("generations", 0),
        "generations_se": data.get("generations_se", 0),
        "years_bp": data.get("years_bp", 0),
        "years_bp_se": data.get("years_bp_se", 0),
        "r_squared": data.get("r_squared", 0),
        "amplitude": data.get("amplitude", 0),
    }
