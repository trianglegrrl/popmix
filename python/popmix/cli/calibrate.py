"""Calibration by simulation: a design's rejection rate over simulated replicates (TASKS.md H7)."""

from __future__ import annotations

import textwrap

import click

import popmix as pm

from ._common import (
    load_f2_blocks,
    parse_population_list,
    pseudohaploid_option,
    sample_selection_options,
    snp_filter_options,
    write_table,
)
from ._errors import EXPECTED_ERRORS, fail

READING = ("If the simulations make the model true, this is the design's false-positive rate, which should be "
           "near alpha; if they make it false, it is its power.")


def _replicate_options(command):
    for option in reversed([
        click.argument("replicates", nargs=-1, required=True),
        click.option("--alpha", type=click.FloatRange(0, 1, min_open=True, max_open=True), default=0.05,
                     help="Significance level: a replicate rejects at p at or below it (default: 0.05)"),
        click.option("--format", "output_format", type=click.Choice(["text", "tsv", "json"]), default="text",
                     help="text: the summary; tsv, json: one row per replicate"),
        click.option("--output", "-o", default=None, help="Also write one row per replicate as TSV here"),
        click.option("--block-size", type=float, default=0.05, help="Block size in Morgans (default: 0.05)"),
        click.option("--maxmiss", type=click.FloatRange(0, 1), default=0.0,
                     help="Drop SNPs missing in more than this fraction of populations (default: 0)"),
        click.option("--all-chromosomes", is_flag=True, default=False,
                     help="Use all chromosomes, not only autosomes 1-22"),
        pseudohaploid_option,
        *snp_filter_options,
        *sample_selection_options,
    ]):
        command = option(command)
    return command


def _f2_list(replicates: tuple, pops: list, block_size: float, maxmiss: float, all_chromosomes: bool) -> list:
    return [load_f2_blocks(r, None, pops, block_size, maxmiss, all_chromosomes) for r in replicates]


def _report(df, replicates: tuple, title: str, output_format: str, output: str | None) -> None:
    df["replicate"] = list(replicates)
    if output:
        write_table(df, "tsv", output)
    if output_format != "text":
        write_table(df, output_format, None)
        return
    a = df.attrs
    click.echo(title)
    click.echo(f"{a['n']} replicates, {a['rejected']} rejected at alpha {a['alpha']:g}: rate {a['rate']:.3f} "
               f"(95% CI {a['ci_low']:.3f}-{a['ci_high']:.3f})")
    click.echo(textwrap.fill(READING, width=79))


@click.command("calibrate-qpadm")
@click.option("--target", required=True, help="Target population")
@click.option("--left", required=True, callback=parse_population_list, help="Source populations (comma-separated)")
@click.option("--right", required=True, callback=parse_population_list, help="Right populations (comma-separated)")
@click.option("--constrained/--unconstrained", default=False, help="Constrain weights to [0, 1]")
@_replicate_options
def calibrate_qpadm(target: str, left: list, right: list, constrained: bool, replicates: tuple, alpha: float,
                    output_format: str, output: str | None, block_size: float, maxmiss: float,
                    all_chromosomes: bool):
    """How often one qpAdm model is rejected over simulated REPLICATES.

    REPLICATES are genotype prefixes or VCF files, one per simulated
    dataset (make them with msprime or another simulator; popmix only
    reads them). Label VCF samples with --relabel.
    """
    try:
        f2s = _f2_list(replicates, [target, *left, *right], block_size, maxmiss, all_chromosomes)
        df = pm.calibrate_qpadm(f2s, left, right, target, alpha=alpha, constrained=constrained)
        _report(df, replicates, f"qpAdm: {target} = mix of {', '.join(left)}; right {', '.join(right)}",
                output_format, output)
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command("calibrate-f4")
@click.option("--pop1", required=True, help="Population 1")
@click.option("--pop2", required=True, help="Population 2")
@click.option("--pop3", required=True, help="Population 3")
@click.option("--pop4", required=True, help="Population 4")
@_replicate_options
def calibrate_f4(pop1: str, pop2: str, pop3: str, pop4: str, replicates: tuple, alpha: float, output_format: str,
                 output: str | None, block_size: float, maxmiss: float, all_chromosomes: bool):
    """How often f4(pop1, pop2; pop3, pop4) is significantly non-zero over simulated REPLICATES.

    REPLICATES are genotype prefixes or VCF files, one per simulated
    dataset; label VCF samples with --relabel.
    """
    try:
        f2s = _f2_list(replicates, [pop1, pop2, pop3, pop4], block_size, maxmiss, all_chromosomes)
        df = pm.calibrate_f4(f2s, pop1, pop2, pop3, pop4, alpha=alpha)
        _report(df, replicates, f"f4({pop1}, {pop2}; {pop3}, {pop4})", output_format, output)
    except EXPECTED_ERRORS as e:
        fail(e)
