"""qpAdm and qpWave commands."""

from __future__ import annotations

import json
import math
import sys

import click

import popmix as pm

from ._common import (
    _json_number,
    _p_text,
    f2_source_options,
    load_f2_blocks,
    parse_population_list,
)
from ._errors import EXPECTED_ERRORS, fail


@click.command()
@click.option("--target", required=True, help="Target population to model")
@click.option("--left", required=True, callback=parse_population_list, help="Source populations (comma-separated)")
@click.option("--right", required=True, callback=parse_population_list, help="Outgroup populations (comma-separated)")
@click.option("--constrained/--unconstrained", default=False, help="Constrain weights to [0, 1]")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpadm(target: str, left: list, right: list, constrained: bool, input_prefix: str | None,
          f2_dir: str | None, block_size: float, maxmiss: float, all_chromosomes: bool,
          output_format: str):
    """Estimate admixture weights using qpAdm.

    Models TARGET as a mixture of LEFT (source) populations, using RIGHT
    populations as outgroups to constrain the model.
    """
    if len(left) < 1:
        click.echo("Error: Need at least 1 source population in --left", err=True)
        sys.exit(1)
    if len(right) < 2:
        click.echo("Error: Need at least 2 outgroup populations in --right", err=True)
        sys.exit(1)

    all_pops = [target] + list(left) + list(right)
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, all_pops, block_size, maxmiss, all_chromosomes)

        result = pm.qpadm(f2_blocks, left, right, target, constrained=constrained)

        if output_format == "json":
            output = {
                "target": result.target,
                "sources": result.sources,
                "weights": dict(zip(result.sources, result.weights)),
                "se": dict(zip(result.sources, result.weight_se)),
                "chisq": result.chisq,
                "dof": result.dof,
                "p": _json_number(result.p),
                "feasible": result.feasible,
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("target\tsource\tweight\tchisq\tdof\tp")
            for source, weight in zip(result.sources, result.weights):
                click.echo(f"{target}\t{source}\t{weight:.4f}\t{result.chisq:.4f}\t{result.dof}\t{_p_text(result.p)}")
        else:
            click.echo(f"qpAdm: {target} = mix of {', '.join(left)}")
            click.echo(f"Outgroups: {', '.join(right)}")
            click.echo()
            click.echo("Admixture weights:")
            for source, weight in zip(result.sources, result.weights):
                click.echo(f"  {source}: {weight:.4f} ({weight*100:.1f}%)")
            click.echo()
            click.echo(f"Model fit: chisq={result.chisq:.4f}, dof={result.dof}, p={_p_text(result.p)}")
            click.echo(f"Feasible: {result.feasible}")
            if math.isnan(result.p):
                click.echo("=> Model not tested: it has 0 degrees of freedom (as many right populations as"
                           " left), so the data cannot reject it")
            elif result.is_good_fit():
                click.echo("=> Good model fit (p > 0.05)")
            else:
                click.echo("=> Poor model fit (p <= 0.05)")

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--left", required=True, callback=parse_population_list, help="Left populations (comma-separated)")
@click.option("--right", required=True, callback=parse_population_list, help="Right populations (comma-separated)")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpwave(left: list, right: list, input_prefix: str | None, f2_dir: str | None,
           block_size: float, maxmiss: float, all_chromosomes: bool, output_format: str):
    """Test admixture wave count using qpWave.

    Tests whether LEFT populations can be modeled as mixtures of common
    ancestral sources, using RIGHT populations as outgroups.
    """
    if len(left) < 2:
        click.echo("Error: Need at least 2 populations in --left", err=True)
        sys.exit(1)
    if len(right) < 2:
        click.echo("Error: Need at least 2 populations in --right", err=True)
        sys.exit(1)

    all_pops = list(left) + list(right)
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, all_pops, block_size, maxmiss, all_chromosomes)

        result = pm.qpwave(f2_blocks, left, right)

        if output_format == "json":
            output = {
                "left": result.left,
                "right": result.right,
                "estimated_rank": result.estimated_rank,
                "rank_tests": [
                    {"rank": t.rank, "chisq": t.chisq, "dof": t.dof, "p": _json_number(t.p)}
                    for t in result.rank_tests
                ],
                "is_clade": result.is_clade(),
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("rank\tchisq\tdof\tp")
            for t in result.rank_tests:
                click.echo(f"{t.rank}\t{t.chisq:.4f}\t{t.dof}\t{_p_text(t.p, 6)}")
        else:
            click.echo(f"qpWave: {', '.join(left)}")
            click.echo(f"Outgroups: {', '.join(right)}")
            click.echo()
            click.echo("Rank tests:")
            click.echo("  Rank   chisq    dof          p")
            for t in result.rank_tests:
                sig = "*" if t.p < 0.05 else " "
                click.echo(f"  {t.rank:4d}  {t.chisq:8.4f}  {t.dof:4d}  {_p_text(t.p, 6):>8} {sig}")
            click.echo()
            click.echo(f"Estimated rank: {result.estimated_rank}")
            if result.is_clade():
                click.echo("=> Populations form a clade (rank 0 not rejected)")
            else:
                click.echo(f"=> {result.estimated_rank + 1} independent ancestry streams")

    except EXPECTED_ERRORS as e:
        fail(e)
