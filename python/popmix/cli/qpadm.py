"""qpAdm and qpWave commands."""

from __future__ import annotations

import dataclasses
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
@click.option("--compat", type=click.Choice(["admixtools2-2.0.10"]), default=None,
              help="Reproduce admixtools2's popdrop nested path, which compares siblings (see DEVIATIONS.md)")
@click.option("--fudge", type=float, default=1e-4,
              help="Regularisation of the f4 covariance, added to its diagonal as a fraction of its trace "
                   "(admixtools2 fudge; default: 1e-4)")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpadm(target: str, left: list, right: list, constrained: bool, compat: str | None, fudge: float,
          input_prefix: str | None,
          f2_dir: str | None, block_size: float, maxmiss: float, all_chromosomes: bool,
          output_format: str):
    """Estimate admixture weights using qpAdm.

    Models TARGET as a mixture of LEFT (source) populations, using RIGHT
    populations as outgroups to constrain the model. Also fits every subset
    of the sources (admixtools2's popdrop table: pattern 1 marks a source
    left out; best marks the nested path).
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

        result = pm.qpadm(f2_blocks, left, right, target, constrained=constrained, fudge=fudge, compat=compat)

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
                "popdrop": [_popdrop_json(row) for row in result.popdrop],
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("target\tsource\tweight\tse\tchisq\tdof\tp")
            for source, weight, se in zip(result.sources, result.weights, result.weight_se):
                click.echo(f"{target}\t{source}\t{weight:.4f}\t{se:.4f}\t{result.chisq:.4f}\t{result.dof}"
                           f"\t{_p_text(result.p)}")
        else:
            click.echo(f"qpAdm: {target} = mix of {', '.join(left)}")
            click.echo(f"Outgroups: {', '.join(right)}")
            click.echo()
            click.echo("Admixture weights:")
            for source, weight, se in zip(result.sources, result.weights, result.weight_se):
                click.echo(f"  {source}: {weight:.4f} ({weight*100:.1f}%), SE {se:.4f}")
            click.echo()
            click.echo(f"Model fit: chisq={result.chisq:.4f}, dof={result.dof}, p={_p_text(result.p)}")
            click.echo(f"Feasible: {result.feasible}")
            if math.isnan(result.p):
                click.echo("=> Model not tested: it has 0 degrees of freedom (as many right populations as"
                           " left), so the data cannot reject it")
            elif result.p > 0.05:  # G1 replaces this verdict
                click.echo("=> Good model fit (p > 0.05)")
            else:
                click.echo("=> Poor model fit (p <= 0.05)")
            _echo_popdrop(result.popdrop)

    except EXPECTED_ERRORS as e:
        fail(e)


def _popdrop_json(row: pm.PopDropRow) -> dict:
    out = dataclasses.asdict(row)
    for key in ("chisq", "p", "chisqdiff", "p_nested"):
        out[key] = _json_number(out[key]) if out[key] is not None else None
    out["weights"] = [None if w is None else _json_number(w) for w in row.weights]
    return out


def _echo_popdrop(rows: tuple[pm.PopDropRow, ...]) -> None:
    """The popdrop table: each subset of the sources, as admixtools2 prints it."""
    click.echo()
    click.echo("Models (popdrop; pattern 1 = source left out):")
    click.echo(f"  {'pattern':<8} {'wt':>3} {'dof':>4} {'chisq':>10} {'p':>9} {'f4rank':>6} {'feasible':>8} "
               f"{'best':>5} weights")
    for r in rows:
        best = "" if r.best is None else str(r.best)
        weights = ", ".join("NA" if w is None else f"{w:.4f}" for w in r.weights)
        click.echo(f"  {r.pattern:<8} {r.wt:>3} {r.dof:>4} {r.chisq:>10.4f} {_p_text(r.p):>9} {r.f4rank:>6} "
                   f"{str(r.feasible):>8} {best:>5} {weights}")


@click.command()
@click.option("--left", required=True, callback=parse_population_list, help="Left populations (comma-separated)")
@click.option("--right", required=True, callback=parse_population_list, help="Right populations (comma-separated)")
@click.option("--fudge", type=float, default=1e-4,
              help="Regularisation of the f4 covariance, added to its diagonal as a fraction of its trace "
                   "(admixtools2 fudge; default: 1e-4)")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpwave(left: list, right: list, fudge: float, input_prefix: str | None, f2_dir: str | None,
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

        result = pm.qpwave(f2_blocks, left, right, fudge=fudge)

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
