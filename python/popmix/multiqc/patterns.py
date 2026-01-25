"""File pattern matching for MultiQC popmix module.

Defines search patterns for finding popmix output files.
"""

import os
from pathlib import Path
from typing import Dict, List, Any


# File patterns for popmix output files
POPMIX_PATTERNS = [
    "*_popmix_stats.json",
    "*_qpadm.json",
    "*_qpwave.json",
    "*_dates.json",
    "*_f2.json",
    "*_f3.json",
    "*_f4.json",
    "*_fst.json",
]


def find_popmix_files(directory: str) -> List[str]:
    """Find all popmix output files in a directory.
    
    Parameters
    ----------
    directory : str
        Directory to search
    
    Returns
    -------
    List[str]
        List of file paths
    """
    found_files = []
    dir_path = Path(directory)
    
    for pattern in POPMIX_PATTERNS:
        found_files.extend(str(f) for f in dir_path.glob(pattern))
    
    # Also search recursively
    for pattern in POPMIX_PATTERNS:
        found_files.extend(str(f) for f in dir_path.rglob(pattern))
    
    # Remove duplicates while preserving order
    seen = set()
    unique_files = []
    for f in found_files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)
    
    return unique_files


def get_search_patterns() -> Dict[str, Dict[str, Any]]:
    """Get MultiQC search patterns configuration.
    
    Returns
    -------
    Dict[str, Dict[str, Any]]
        Search patterns dictionary for MultiQC
    """
    return {
        "popmix": {
            "fn": "*_popmix_stats.json",
        },
        "popmix/qpadm": {
            "fn": "*_qpadm.json",
            "contents": '"target"',
        },
        "popmix/qpwave": {
            "fn": "*_qpwave.json",
            "contents": '"rank_tests"',
        },
        "popmix/dates": {
            "fn": "*_dates.json",
            "contents": '"generations"',
        },
        "popmix/f2": {
            "fn": "*_f2.json",
        },
        "popmix/f3": {
            "fn": "*_f3.json",
        },
        "popmix/f4": {
            "fn": "*_f4.json",
        },
        "popmix/fst": {
            "fn": "*_fst.json",
        },
    }


def get_sample_name(filepath: str) -> str:
    """Extract sample name from file path.
    
    Parameters
    ----------
    filepath : str
        Path to popmix output file
    
    Returns
    -------
    str
        Sample name
    """
    filename = Path(filepath).stem
    
    # Remove common suffixes
    suffixes = [
        "_popmix_stats",
        "_qpadm",
        "_qpwave",
        "_dates",
        "_f2",
        "_f3",
        "_f4",
        "_fst",
    ]
    
    for suffix in suffixes:
        if filename.endswith(suffix):
            return filename[:-len(suffix)]
    
    return filename
