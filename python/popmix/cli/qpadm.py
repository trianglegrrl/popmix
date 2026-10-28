"""qpAdm and qpWave commands."""

from __future__ import annotations

import json
import sys
import textwrap

import click

import popmix as pm

from . import _qpadm_output as out
from ._common import (
    _json_number,
    _p_text,
    _pseudohaploid,
    f2_source_options,
    load_f2_blocks,
    parse_population_list,
    sample_selection,
    snp_filters,
    write_columns,
)
from ._errors import EXPECTED_ERRORS, fail
from ._population_checks import validate_populations


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
@click.option("--alpha", type=click.FloatRange(0, 1, min_open=True, max_open=True), default=0.05,
              help="Significance level: p at or below it rejects the model (default: 0.05)")
@f2_source_options
@click.option("--f4-snps", type=click.Choice(["f2", "model", "each"]), default="f2",
              help="With --input: where each f4 comes from. f2 (default): f2 blocks on the SNPs every population has data for, polymorphic ones, as admixtools2 qpadm(f2_from_geno(...)); model: from the genotypes on the SNPs every population has data for, as admixtools2 qpadm(prefix); each: each f4 on the SNPs its four populations have data for, as "
                   "qpadm(prefix, allsnps = TRUE)")
@click.option("--output-prefix", "-o", default=None,
              help="Also write PREFIX.weights.tsv, .rankdrop.tsv, .popdrop.tsv and .json")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpadm(target: str, left: list, right: list, constrained: bool, compat: str | None, fudge: float,
          alpha: float, input_prefix: str | None,
          f2_dir: str | None, block_size: float, maxmiss: float, all_chromosomes: bool, f4_snps: str,
          output_prefix: str | None, output_format: str):
    """Estimate admixture weights using qpAdm.

    Models TARGET as a mixture of LEFT (source) populations, using RIGHT
    populations as outgroups to constrain the model. Also gives the rank
    tests of the target and sources (admixtools2's rankdrop) and fits every
    subset of the sources (its popdrop table: pattern 1 marks a source left
    out; best marks the nested path).
    """
    if len(left) < 1:
        click.echo("Error: Need at least 1 source population in --left", err=True)
        sys.exit(1)
    if len(right) < 2:
        click.echo("Error: Need at least 2 outgroup populations in --right", err=True)
        sys.exit(1)

    all_pops = [target] + list(left) + list(right)
    try:
        data, kw = _f4_source(f4_snps, input_prefix, f2_dir, all_pops, block_size, maxmiss, all_chromosomes)
        result = pm.qpadm(data, left, right, target, constrained=constrained, fudge=fudge, compat=compat, **kw)
        assessment = result.assess(alpha)

        if output_prefix:
            write_columns(out.weight_columns(result), f"{output_prefix}.weights.tsv")
            write_columns(out.rankdrop_columns(result), f"{output_prefix}.rankdrop.tsv")
            write_columns(out.popdrop_columns(result), f"{output_prefix}.popdrop.tsv")
            with open(f"{output_prefix}.json", "w") as f:
                json.dump(out.as_json(result, assessment), f, indent=2)
                f.write("\n")
        if output_format == "json":
            click.echo(json.dumps(out.as_json(result, assessment), indent=2))
        elif output_format == "tsv":
            write_columns(out.weight_columns(result), None)
        else:
            _echo_text(result, right, assessment)

    except EXPECTED_ERRORS as e:
        fail(e)


def _f4_source(f4_snps: str, input_prefix: str | None, f2_dir: str | None, pops: list, block_size: float,
               maxmiss: float, all_chromosomes: bool) -> tuple:
    """The data qpAdm or qpWave fits, and its keywords: f2 blocks, or with
    --f4-snps model or each, the genotype prefix (TASKS.md H11)."""
    if f4_snps == "f2":
        return load_f2_blocks(input_prefix, f2_dir, pops, block_size, maxmiss, all_chromosomes), {}
    if input_prefix is None:
        raise click.UsageError(f"--f4-snps {f4_snps} computes f4 from genotypes: give --input, not --f2-dir")
    given = [flag for flag, on in [("--maxmiss", maxmiss != 0), ("SNP filters", bool(snp_filters())),
                                   ("sample selection", bool(sample_selection()))] if on]
    if given:
        raise click.UsageError(f"{', '.join(given)} apply to f2 blocks, not to --f4-snps {f4_snps}")
    validate_populations(input_prefix, pops)
    return input_prefix, {"allsnps": f4_snps == "each", "block_size": block_size,
                          "auto_only": not all_chromosomes, "pseudohaploid": _pseudohaploid()}


def _echo_text(result: pm.QpAdmResult, right: list, assessment: pm.QpAdmAssessment) -> None:
    click.echo(f"qpAdm: {result.target} = mix of {', '.join(result.sources)}")
    click.echo(f"Outgroups: {', '.join(right)}")
    click.echo(f"{result.n_blocks} jackknife blocks, {result.n_snps} SNPs")
    click.echo()
    click.echo("Admixture weights:")
    width = max(len(s) for s in result.sources)
    click.echo(f"  {'source':<{width}} {'weight':>8} {'SE':>9} {'z':>8}")
    for source, weight, se, z in zip(result.sources, result.weights, result.weight_se, result.weight_z):
        click.echo(f"  {source:<{width}} {weight:>8.4f} {se:>9.3g} {z:>8.1f}")
    click.echo()
    click.echo(f"Model fit: chisq={result.chisq:.4f}, dof={result.dof}, p={out.p_text(result.p)}")
    click.echo(f"Feasible: {result.feasible}")
    _echo_assessment(assessment)
    click.echo()
    click.echo("Rank tests (rankdrop):")
    click.echo(f"  {'f4rank':>6} {'dof':>4} {'chisq':>12} {'p':>10}")
    for t in result.rankdrop:
        click.echo(f"  {t.rank:>6} {t.dof:>4} {t.chisq:>12.4f} {out.p_text(t.p):>10}")
    _echo_popdrop(result.popdrop)


_VERDICTS = {
    "not tested": "=> Not tested: 0 degrees of freedom (as many right populations as left), so the data cannot "
                  "reject the model.",
    "rejected": "=> Rejected at alpha {alpha:g}: the sources do not explain the target's f4 relations to the "
                "right populations.",
    "infeasible": "=> Rejected: p is above alpha {alpha:g}, but a weight is outside [0, 1].",
    "not rejected": "=> Not rejected at alpha {alpha:g}, with every weight in [0, 1]. That does not make it the "
                    "model: others may fit too.",
}


def _echo_assessment(assessment: pm.QpAdmAssessment) -> None:
    """The verdict from p, feasibility, SEs and popdrop, never from p alone (G1)."""
    click.echo(textwrap.fill(_VERDICTS[assessment.status].format(alpha=assessment.alpha), width=79,
                             subsequent_indent="   "))
    for note in assessment.notes:
        click.echo(textwrap.fill(note, width=79, initial_indent="   - ", subsequent_indent="     "))


def _echo_popdrop(rows: tuple[pm.PopDropRow, ...]) -> None:
    """The popdrop table: each subset of the sources, as admixtools2 prints it."""
    click.echo()
    click.echo("Models (popdrop; pattern 1 = source left out):")
    click.echo(f"  {'pattern':<8} {'wt':>3} {'dof':>4} {'chisq':>10} {'p':>9} {'f4rank':>6} {'feasible':>8} "
               f"{'best':>5} weights")
    for r in rows:
        best = "" if r.best is None else str(r.best)
        weights = ", ".join("NA" if w is None else f"{w:.4f}" for w in r.weights)
        click.echo(f"  {r.pattern:<8} {r.wt:>3} {r.dof:>4} {r.chisq:>10.4f} {out.p_text(r.p):>9} {r.f4rank:>6} "
                   f"{str(r.feasible):>8} {best:>5} {weights}")


def _rank_reading(rank: int, max_rank: int) -> str:
    """What qpWave's estimated rank shows: a lower bound, never a clade (G2)."""
    more = f"other right populations or more data can show more streams (the rank can be at most {max_rank} here)."
    if rank == 0:
        text = ("=> Rank 0 not rejected: these right populations do not distinguish the left ones. "
                f"That does not show that they form a clade: {more}")
    else:
        text = (f"=> At least {rank + 1} streams of ancestry, as seen from these right populations: a lower "
                f"bound, as {more}")
    return textwrap.fill(text, width=79, subsequent_indent="   ")


@click.command()
@click.option("--left", required=True, callback=parse_population_list, help="Left populations (comma-separated)")
@click.option("--right", required=True, callback=parse_population_list, help="Right populations (comma-separated)")
@click.option("--fudge", type=float, default=1e-4,
              help="Regularisation of the f4 covariance, added to its diagonal as a fraction of its trace "
                   "(admixtools2 fudge; default: 1e-4)")
@click.option("--alpha", type=click.FloatRange(0, 1, min_open=True, max_open=True), default=0.05,
              help="Significance level of the rank tests (default: 0.05)")
@f2_source_options
@click.option("--f4-snps", type=click.Choice(["f2", "model"]), default="f2",
              help="With --input: where each f4 comes from. f2 (default): f2 blocks on the SNPs every population has data for, polymorphic ones, as admixtools2 qpadm(f2_from_geno(...)); model: from the genotypes on the SNPs every population has data for, as admixtools2 qpadm(prefix)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpwave(left: list, right: list, fudge: float, alpha: float, input_prefix: str | None, f2_dir: str | None,
           block_size: float, maxmiss: float, all_chromosomes: bool, f4_snps: str, output_format: str):
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
        data, kw = _f4_source(f4_snps, input_prefix, f2_dir, all_pops, block_size, maxmiss, all_chromosomes)
        kw.pop("allsnps", None)
        result = pm.qpwave(data, left, right, fudge=fudge, alpha=alpha, **kw)

        if output_format == "json":
            output = {
                "left": result.left,
                "right": result.right,
                "alpha": result.alpha,
                "estimated_rank": result.estimated_rank,
                "rank_tests": [
                    {"rank": t.rank, "chisq": t.chisq, "dof": t.dof, "p": _json_number(t.p)}
                    for t in result.rank_tests
                ],
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
            click.echo(f"Rank tests (* rejected at alpha {alpha:g}):")
            click.echo("  Rank   chisq    dof          p")
            for t in result.rank_tests:
                sig = "*" if t.p <= alpha else " "
                click.echo(f"  {t.rank:4d}  {t.chisq:8.4f}  {t.dof:4d}  {_p_text(t.p, 6):>8} {sig}")
            click.echo()
            click.echo(f"Estimated rank: {result.estimated_rank} (the lowest rank not rejected)")
            click.echo(_rank_reading(result.estimated_rank, min(len(left), len(right)) - 1))

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command("qpadm-rotate")
@click.option("--target", required=True, help="Target population to model")
@click.option("--leftright", required=True, callback=parse_population_list,
              help="Populations distributed between sources and right (comma-separated)")
@click.option("--rightfix", callback=parse_population_list, default="",
              help="Populations on the right in every model (comma-separated)")
@click.option("--constrained/--unconstrained", default=False, help="Constrain weights to [0, 1]")
@click.option("--fudge", type=float, default=1e-4,
              help="Regularisation of the f4 covariance (admixtools2 fudge; default: 1e-4)")
@click.option("--alpha", type=click.FloatRange(0, 1, min_open=True, max_open=True), default=0.05,
              help="Significance level for each model's status (default: 0.05)")
@f2_source_options
@click.option("--output", "-o", default=None, help="Write the table here instead of stdout")
def qpadm_rotate(target: str, leftright: list, rightfix: list, constrained: bool, fudge: float, alpha: float,
                 input_prefix: str | None, f2_dir: str | None, block_size: float, maxmiss: float,
                 all_chromosomes: bool, output: str | None):
    """Fit every qpAdm model of a rotation (admixtools2 qpadm_rotate).

    Each subset of LEFTRIGHT in turn is the sources, and the rest of it,
    after RIGHTFIX, the right populations. Writes one TSV row per model:
    its sources and right populations, chisq, dof, p, feasibility, status
    at --alpha, and the weights and SEs in source order. The p-values are
    not corrected for testing every model; the number tested is warned.
    """
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, [target, *leftright, *rightfix], block_size, maxmiss,
                                   all_chromosomes)
        fits = pm.qpadm_rotate(f2_blocks, leftright, target, rightfix or None, full_results=True, fudge=fudge,
                               constrained=constrained)
        def joined(xs) -> str:
            return ",".join(f"{x:.12g}" for x in xs)

        write_columns({
            "model": list(range(1, len(fits) + 1)), "target": [target] * len(fits),
            "left": [",".join(r.sources) for r in fits], "right": [",".join(r.right) for r in fits],
            "n_sources": [len(r.sources) for r in fits], "chisq": [r.chisq for r in fits],
            "dof": [r.dof for r in fits], "p": [r.p for r in fits], "feasible": [r.feasible for r in fits],
            "status": [r.assess(alpha).status for r in fits], "weights": [joined(r.weights) for r in fits],
            "se": [joined(r.weight_se) for r in fits], "n_blocks": [r.n_blocks for r in fits],
            "n_snps": [r.n_snps for r in fits],
        }, output)
    except EXPECTED_ERRORS as e:
        fail(e)

