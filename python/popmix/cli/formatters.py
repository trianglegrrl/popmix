"""Output formatters for popmix results.

This module provides machine-readable output formats for pipeline integration:
- JSON: Structured, nested data
- YAML: Human-readable configuration format
- TSV: Tab-separated for spreadsheets
- NDJSON: Newline-delimited JSON for streaming
- Parquet: High-performance columnar format
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# Import YAML conditionally (it's an optional dependency)
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

# Import pandas/pyarrow conditionally for Parquet support
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False


# =============================================================================
# Confidence Interval Calculation
# =============================================================================

def calculate_confidence_interval(
    estimate: float,
    std_err: float,
    confidence: float = 0.95,
) -> tuple:
    """Calculate confidence interval for an estimate.
    
    Parameters
    ----------
    estimate : float
        Point estimate
    std_err : float
        Standard error of the estimate
    confidence : float
        Confidence level (default 0.95)
    
    Returns
    -------
    tuple
        (lower, upper) bounds
    """
    z_scores = {0.90: 1.645, 0.95: 1.96, 0.99: 2.576}
    z = z_scores.get(confidence, 1.96)
    margin = z * std_err
    return (estimate - margin, estimate + margin)


# =============================================================================
# Result to Dictionary Conversion
# =============================================================================

def result_to_dict(
    result: Any,
    result_type: str,
    include_metadata: bool = False,
    include_ci: bool = True,
) -> Dict[str, Any]:
    """Convert a result object to a dictionary.
    
    Parameters
    ----------
    result : Any
        Result object (F2, F3, F4, qpAdm, DATES, etc.)
    result_type : str
        Type of result: "f2", "f3", "f4", "qpadm", "qpwave", "dates"
    include_metadata : bool
        Whether to include result type metadata
    include_ci : bool
        Whether to include confidence intervals
    
    Returns
    -------
    Dict[str, Any]
        Dictionary representation
    """
    data = {}
    
    if include_metadata:
        data["type"] = result_type
        data["result_type"] = result_type
    
    if result_type == "f2":
        data.update({
            "pop1": result.pop1,
            "pop2": result.pop2,
            "estimate": result.estimate,
            "std_err": result.std_err,
            "z_score": result.z_score,
            "n_snps": result.n_snps,
        })
        
        if include_ci and hasattr(result, 'std_err'):
            lower, upper = calculate_confidence_interval(result.estimate, result.std_err)
            data["ci_lower"] = lower
            data["ci_upper"] = upper
            data["confidence_interval"] = [lower, upper]
    
    elif result_type == "f3":
        data.update({
            "target": result.target,
            "source1": result.source1,
            "source2": result.source2,
            "estimate": result.estimate,
            "std_err": result.std_err,
            "z_score": result.z_score,
            "n_blocks": result.n_blocks,
        })
        
        if include_ci and hasattr(result, 'std_err'):
            lower, upper = calculate_confidence_interval(result.estimate, result.std_err)
            data["ci_lower"] = lower
            data["ci_upper"] = upper
    
    elif result_type == "f4":
        data.update({
            "pop1": result.pop1,
            "pop2": result.pop2,
            "pop3": result.pop3,
            "pop4": result.pop4,
            "estimate": result.estimate,
            "std_err": result.std_err,
            "z_score": result.z_score,
            "n_blocks": result.n_blocks,
        })
        
        if include_ci and hasattr(result, 'std_err'):
            lower, upper = calculate_confidence_interval(result.estimate, result.std_err)
            data["ci_lower"] = lower
            data["ci_upper"] = upper
    
    elif result_type == "qpadm":
        # Build weights dictionary
        weights_dict = {}
        weights_list = []
        
        for i, source in enumerate(result.sources):
            weight = result.weights[i]
            se = result.std_errs[i] if hasattr(result, 'std_errs') and result.std_errs else None
            
            weight_data = {"source": source, "weight": weight}
            if se is not None:
                weight_data["std_err"] = se
                lower, upper = calculate_confidence_interval(weight, se)
                weight_data["ci_lower"] = lower
                weight_data["ci_upper"] = upper
            
            weights_dict[source] = weight
            weights_list.append(weight_data)
        
        # Check feasibility
        feasible = all(0 <= w <= 1 for w in result.weights)
        good_fit = result.p_value > 0.05
        
        data.update({
            "target": result.target,
            "sources": result.sources,
            "weights": weights_dict,
            "weights_detail": weights_list,
            "chi_sq": result.chi_sq,
            "df": result.df,
            "p_value": result.p_value,
            "model_fit": {
                "chi_sq": result.chi_sq,
                "df": result.df,
                "p_value": result.p_value,
                "good_fit": good_fit,
            },
            "feasible": feasible,
        })
    
    elif result_type == "qpwave":
        rank_tests = []
        for test in result.rank_tests:
            rank_tests.append({
                "rank": test.rank,
                "chi_sq": test.chi_sq,
                "df": test.df,
                "p_value": test.p_value,
                "rejected": test.p_value < 0.05,
            })
        
        data.update({
            "left_pops": result.left_pops,
            "right_pops": result.right_pops,
            "estimated_rank": result.estimated_rank,
            "rank_tests": rank_tests,
            "is_clade": result.is_clade() if hasattr(result, 'is_clade') else None,
        })
    
    elif result_type == "dates":
        data.update({
            "test_pop": result.test_pop,
            "source1": result.source1,
            "source2": result.source2,
            "generations": result.generations,
            "generations_se": result.generations_se,
            "years_bp": result.years_bp,
            "years_bp_se": result.years_bp_se,
            "amplitude": result.amplitude,
            "r_squared": result.r_squared,
            "n_snps": result.n_snps,
            "timing": {
                "generations": result.generations,
                "generations_se": result.generations_se,
                "years_bp": result.years_bp,
                "years_bp_se": result.years_bp_se,
            },
            "fit": {
                "amplitude": result.amplitude,
                "r_squared": result.r_squared,
            },
        })
    
    else:
        # Generic serialization for unknown types
        for attr in dir(result):
            if not attr.startswith('_') and not callable(getattr(result, attr)):
                data[attr] = getattr(result, attr)
    
    return data


# =============================================================================
# JSON Formatting
# =============================================================================

def format_as_json(
    result: Any,
    result_type: str,
    indent: int = 2,
    include_metadata: bool = False,
) -> str:
    """Format result as JSON string.
    
    Parameters
    ----------
    result : Any
        Result object
    result_type : str
        Type of result
    indent : int
        JSON indentation level
    include_metadata : bool
        Whether to include metadata
    
    Returns
    -------
    str
        JSON string
    """
    data = result_to_dict(result, result_type, include_metadata=include_metadata)
    return json.dumps(data, indent=indent)


# =============================================================================
# YAML Formatting
# =============================================================================

def format_as_yaml(
    result: Any,
    result_type: str,
    include_metadata: bool = False,
) -> str:
    """Format result as YAML string.
    
    Parameters
    ----------
    result : Any
        Result object
    result_type : str
        Type of result
    include_metadata : bool
        Whether to include metadata
    
    Returns
    -------
    str
        YAML string
    """
    if not YAML_AVAILABLE:
        raise ImportError("PyYAML is required for YAML output. Install with: pip install pyyaml")
    
    data = result_to_dict(result, result_type, include_metadata=include_metadata)
    return yaml.dump(data, default_flow_style=False, sort_keys=False)


# =============================================================================
# TSV Formatting
# =============================================================================

def format_as_tsv(
    result: Any,
    result_type: str,
) -> str:
    """Format result as TSV (tab-separated values) string.
    
    Parameters
    ----------
    result : Any
        Result object
    result_type : str
        Type of result
    
    Returns
    -------
    str
        TSV string with header
    """
    data = result_to_dict(result, result_type, include_ci=True)
    
    # Flatten nested structures
    flat_data = {}
    for key, value in data.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                flat_data[f"{key}_{sub_key}"] = sub_value
        elif isinstance(value, list) and len(value) > 0 and isinstance(value[0], dict):
            # Skip complex nested lists in TSV
            continue
        elif isinstance(value, list):
            flat_data[key] = ",".join(str(v) for v in value)
        else:
            flat_data[key] = value
    
    # Build TSV
    headers = list(flat_data.keys())
    values = [str(flat_data[h]) for h in headers]
    
    return "\t".join(headers) + "\n" + "\t".join(values)


# =============================================================================
# NDJSON Formatting
# =============================================================================

def format_as_ndjson(
    results: List[Any],
    result_type: str,
) -> str:
    """Format multiple results as NDJSON (newline-delimited JSON).
    
    Parameters
    ----------
    results : List[Any]
        List of result objects
    result_type : str
        Type of results
    
    Returns
    -------
    str
        NDJSON string (one JSON object per line)
    """
    lines = []
    for result in results:
        data = result_to_dict(result, result_type)
        lines.append(json.dumps(data))
    return "\n".join(lines)


# =============================================================================
# Parquet Formatting
# =============================================================================

def write_as_parquet(
    result: Any,
    result_type: str,
    output_path: str,
) -> None:
    """Write result to Parquet file.
    
    Parameters
    ----------
    result : Any
        Result object
    result_type : str
        Type of result
    output_path : str
        Output file path
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas is required for Parquet output. Install with: pip install pandas pyarrow")
    
    data = result_to_dict(result, result_type)
    
    # Flatten for DataFrame
    flat_data = {}
    for key, value in data.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                if not isinstance(sub_value, (dict, list)):
                    flat_data[f"{key}_{sub_key}"] = sub_value
        elif isinstance(value, list) and len(value) > 0 and isinstance(value[0], dict):
            continue
        elif isinstance(value, list):
            flat_data[key] = str(value)
        else:
            flat_data[key] = value
    
    # Create DataFrame with single row
    df = pd.DataFrame([flat_data])
    
    # Write to Parquet
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, index=False)


# =============================================================================
# File Output
# =============================================================================

def write_result_to_file(
    result: Any,
    result_type: str,
    output_path: str,
    format: Optional[str] = None,
) -> None:
    """Write result to file in specified format.
    
    Parameters
    ----------
    result : Any
        Result object
    result_type : str
        Type of result
    output_path : str
        Output file path
    format : str, optional
        Output format. If None, inferred from file extension.
    """
    path = Path(output_path)
    
    # Create parent directories
    path.parent.mkdir(parents=True, exist_ok=True)
    
    # Infer format from extension if not specified
    if format is None:
        ext = path.suffix.lower()
        format_map = {
            ".json": "json",
            ".yaml": "yaml",
            ".yml": "yaml",
            ".tsv": "tsv",
            ".txt": "tsv",
            ".parquet": "parquet",
            ".ndjson": "ndjson",
            ".jsonl": "ndjson",
        }
        format = format_map.get(ext, "json")
    
    # Write based on format
    if format == "json":
        content = format_as_json(result, result_type)
        path.write_text(content)
    
    elif format == "yaml":
        content = format_as_yaml(result, result_type)
        path.write_text(content)
    
    elif format == "tsv":
        content = format_as_tsv(result, result_type)
        path.write_text(content)
    
    elif format == "parquet":
        write_as_parquet(result, result_type, str(path))
    
    elif format == "ndjson":
        # For single result, wrap in list
        content = format_as_ndjson([result], result_type)
        path.write_text(content)
    
    else:
        raise ValueError(f"Unknown format: {format}")


# =============================================================================
# Format Detection
# =============================================================================

def detect_format_from_path(path: str) -> str:
    """Detect output format from file path.
    
    Parameters
    ----------
    path : str
        File path
    
    Returns
    -------
    str
        Format name
    """
    ext = Path(path).suffix.lower()
    format_map = {
        ".json": "json",
        ".yaml": "yaml",
        ".yml": "yaml",
        ".tsv": "tsv",
        ".txt": "tsv",
        ".parquet": "parquet",
        ".ndjson": "ndjson",
        ".jsonl": "ndjson",
    }
    return format_map.get(ext, "json")
