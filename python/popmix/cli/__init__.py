"""
Popmix CLI - Command-line interface for population genetics analysis.

Usage:
    popmix [OPTIONS] COMMAND [ARGS]...

Commands:
    info        Show dataset information
    populations List all populations with sample counts
    samples     Each sample's SNPs called, heterozygosity and ploidy
    extract-f2  Extract F2 blocks from genotype data
    f2          Compute f2 statistic between two populations
    f3          Compute f3 statistic (admixture test)
    f4          Compute the f4 statistic f4(pop1, pop2; pop3, pop4)
    qpadm       Estimate admixture weights
    qpwave      Test admixture wave count
    qpadm-rotate Fit every qpAdm model of a rotation
    convertf    Convert between genotype file formats
    merge       Merge two genotype datasets
    dates       Estimate admixture date using DATES
    qpgraph     Fit an admixture graph
    find-graphs Search for the best-fitting admixture graphs
    qpdstat     f4- or D-statistics for many population combinations
    qp3pop      f3-statistics for many population combinations
    outgroup-f3 Outgroup f3 matrix, with an optional heatmap
    dstat-matrix D-statistic matrix for pairs of populations, with an optional heatmap
    pca         Principal component analysis as smartpca, with projection
    kinship     Kinship for every pair of samples (KING-robust or READ-style)
    sfs         Site-frequency spectrum (dadi projection or admixtools2 joint_sfs)
    admixture   Ancestry proportions as ADMIXTURE (unsupervised, supervised, projection)
    rolloff     ROLLOFF: admixture dates from weighted LD correlation
    install-skill Install the popmix agent skill for Claude Code
"""

from __future__ import annotations

import click

import popmix as pm

from ._common import write_columns, write_table
from ._errors import print_popmix_warnings
from ._provenance import record
from .datasets import convertf, info, install_skill, merge, populations, samples
from .dating import dates, rolloff
from .fstats import dstat_matrix, extract_f2, f2, f3, f4, outgroup_f3, qp3pop, qpdstat
from .graphs import find_graphs, qpgraph
from .qpadm import qpadm, qpadm_rotate, qpwave
from .structure import admixture, kinship, pca, sfs

__all__ = ["cli", "main", "write_columns", "write_table"]


# with no command, the help on stdout and success, and one usage line, on
# any Click (8.2 changed both: help on stderr, exit 2, an optional COMMAND)
@click.group(invoke_without_command=True, subcommand_metavar="COMMAND [ARGS]...")
@click.option("--threads", "-t", type=int, default=4, help="Number of threads to use (default: 4)")
@click.option("--debug", is_flag=True, help="Show the traceback of an error")
@click.version_option(version=pm.__version__, prog_name="popmix")
@click.pass_context
def cli(ctx, threads: int, debug: bool):
    """Popmix: Population genetics F-statistics and admixture analysis toolkit."""
    ctx.ensure_object(dict)
    ctx.obj["debug"] = debug
    if threads < 1:
        raise click.BadParameter("Thread count must be at least 1")
    ctx.obj["threads"] = threads
    pm.init_threads(threads)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


for _command in (
    info,
    populations,
    samples,
    convertf,
    merge,
    install_skill,
    extract_f2,
    f2,
    f3,
    f4,
    qpdstat,
    qp3pop,
    outgroup_f3,
    dstat_matrix,
    qpadm,
    qpadm_rotate,
    qpwave,
    pca,
    kinship,
    sfs,
    admixture,
    dates,
    rolloff,
    qpgraph,
    find_graphs,
):
    cli.add_command(record(_command))


def main():
    """Entry point for the popmix CLI."""
    print_popmix_warnings()
    cli()



if __name__ == "__main__":
    main()
