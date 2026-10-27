"""Checks of the population names a command is given, with suggestions."""

from __future__ import annotations

import click

import popmix as pm

# Failures of the input (bad names or values, missing or malformed files) or
# of the computation (a singular system): reported in one line. Anything else
# is a bug and keeps its traceback.
EXPECTED_ERRORS = (ValueError, OSError, RuntimeError)


def format_population_error(invalid_pop: str, prefix: str) -> str:
    """Format helpful error message for invalid population name."""
    lines = [f"Population '{invalid_pop}' not found in dataset."]

    try:
        suggestions = pm.suggest_populations(prefix, invalid_pop, limit=5)
        if suggestions:
            lines.append("\nDid you mean one of these?")
            for pop, count in suggestions:
                lines.append(f"  - {pop} ({count} samples)")
    except EXPECTED_ERRORS:
        pass  # no suggestions if the dataset cannot be listed

    lines.append(f"\nTip: List all populations with: popmix populations {prefix}")
    lines.append(f'Tip: Search populations with: popmix populations {prefix} --filter "{invalid_pop[:4]}*"')

    return "\n".join(lines)


def validate_populations(prefix: str, populations: list) -> None:
    """Validate that all populations exist in the dataset.

    Raises click.ClickException with helpful error message if any population is missing.
    """
    if not populations:
        return

    try:
        available = {p for p, _ in pm.list_populations(prefix)}
    except EXPECTED_ERRORS as e:
        click.echo(f"Warning: Could not validate populations: {e}", err=True)
        return

    missing = [p for p in populations if p not in available]
    if missing:
        # Format error message for the first missing population
        error_msg = format_population_error(missing[0], prefix)
        if len(missing) > 1:
            error_msg += f"\n\nAlso missing: {', '.join(missing[1:])}"
        raise click.ClickException(error_msg)
