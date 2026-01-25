"""JSON schemas for popmix result validation."""

from pathlib import Path

SCHEMA_DIR = Path(__file__).parent


def get_schema_path(schema_name: str) -> Path:
    """Get path to a schema file.
    
    Parameters
    ----------
    schema_name : str
        Schema name (e.g., "f2_result")
    
    Returns
    -------
    Path
        Path to schema file
    """
    return SCHEMA_DIR / f"{schema_name}.json"


def load_schema(schema_name: str) -> dict:
    """Load a JSON schema.
    
    Parameters
    ----------
    schema_name : str
        Schema name (e.g., "f2_result")
    
    Returns
    -------
    dict
        Schema dictionary
    """
    import json
    
    schema_path = get_schema_path(schema_name)
    if not schema_path.exists():
        raise FileNotFoundError(f"Schema not found: {schema_path}")
    
    with open(schema_path) as f:
        return json.load(f)
