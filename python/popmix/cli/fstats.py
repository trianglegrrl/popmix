"""f-statistic commands: f2 blocks, f2, f3, f4 and their many-population forms."""

from __future__ import annotations

import json
import sys

import click

import popmix as pm
from popmix import _core

from ._common import (
    _json_row,
    _pseudohaploid,
    f2_source_options,
    load_f2_blocks,
    matrix_options,
    parse_population_list,
    parse_snp_range,
    pseudohaploid_option,
    sample_selection,
    sample_selection_options,
    scan_output_options,
    selected_samples,
    snp_filter_options,
    snp_filters,
    stat_source,
    stat_source_options,
    with_options,
    write_matrix,
    write_table,
)
from ._errors import EXPECTED_ERRORS, fail
from ._population_checks import validate_populations


@click.command("extract-f2")
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--output", "-o", "output_dir", required=True, help="Output directory for F2 blocks")
@click.option("--snp-range", callback=parse_snp_range,
              help="Only SNPs start:end, by index (computed outside admixtools2's extract_f2 pipeline; not with --maxmiss or --pseudohaploid)")
@click.option("--block-size", type=float, default=0.05, help="Block size in Morgans (default: 0.05)")
@click.option("--populations", callback=parse_population_list,
              help="Comma-separated list of populations to include (default: all but Ignore)")
@click.option("--maxmiss", type=click.FloatRange(0, 1), default=0.0,
              help="Drop SNPs missing in more than this fraction of populations (default: 0)")
@click.option("--auto-only/--all-chromosomes", default=True, help="Use only autosomal SNPs (default: auto-only)")
@pseudohaploid_option
@with_options(snp_filter_options)
@with_options(sample_selection_options)
def extract_f2(input_prefix: str, output_dir: str, snp_range: tuple[int, int] | None,
               block_size: float, populations: list, maxmiss: float, auto_only: bool):
    """Extract F2 blocks from genotype data, as admixtools2's extract_f2().

    Computes pairwise F2 statistics for all population pairs and saves them
    for use in downstream analyses (f3, f4, qpAdm, etc. with --f2-dir).
    """
    try:
        if not snp_range:
            if not populations and sample_selection():
                labels = _core.select_samples(input_prefix, **sample_selection())[1]
                # most samples first, as list_populations
                populations = sorted((p for p in set(labels) if p != "Ignore"), key=lambda p: (-labels.count(p), p))
            elif not populations:
                populations = [p for p, _ in pm.list_populations(input_prefix) if p != "Ignore"]
            selected = selected_samples(input_prefix, populations)
            click.echo(f"Streaming F2 extraction from {input_prefix}...")
            click.echo(f"Populations: {', '.join(populations)}")
            f2_blocks = pm.extract_f2(
                input_prefix,
                populations,
                auto_only=auto_only,
                block_size=block_size,
                maxmiss=maxmiss,
                pseudohaploid=_pseudohaploid(),
                **snp_filters(),
                **selected,
            )
            click.echo(f"Computed {f2_blocks.n_blocks} blocks for {f2_blocks.n_populations} populations")
        else:
            # A SNP index range needs the genotypes in memory, a path without
            # the missingness filter or a choice of ploidy rule
            for flag, unset in [("--maxmiss", maxmiss == 0), ("--pseudohaploid", _pseudohaploid() == "detect"),
                                ("SNP filters", not snp_filters()), ("sample selection", not sample_selection())]:
                if not unset:
                    raise click.UsageError(f"{flag} does not apply with --snp-range")
            click.echo(f"Reading genotypes from {input_prefix}...")
            geno = pm.read_genotypes(input_prefix, snp_range=snp_range)
            click.echo(f"Loaded {geno.n_snps} SNPs x {geno.n_samples} samples")
            if auto_only:
                other = sum(not (s.chromosome.isdigit() and 1 <= int(s.chromosome) <= 22) for s in geno.snps)
                if other:
                    raise click.UsageError(f"the SNP range holds {other} non-autosomal SNPs: narrow it, or "
                                           "pass --all-chromosomes to use them")
            if populations:
                validate_populations(input_prefix, populations)
            else:
                populations = [p for p in geno.populations if p != "Ignore"]

            click.echo("Computing allele frequencies...")
            af = geno.allele_frequencies(populations=populations)
            click.echo(f"Found {af.n_populations} populations")

            click.echo("Computing F2 blocks...")
            f2_blocks = pm.f2_from_frequencies(af, block_size=block_size)
            click.echo(f"Computed {f2_blocks.n_blocks} blocks")

        # Save F2 blocks (creates directory, metadata.json, and .npy files)
        click.echo(f"Saving F2 blocks to {output_dir}...")
        f2_blocks.save(output_dir, source=input_prefix)
        click.echo(f"F2 blocks saved to {output_dir}")

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--pop1", required=True, help="First population")
@click.option("--pop2", required=True, help="Second population")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def f2(pop1: str, pop2: str, input_prefix: str | None, f2_dir: str | None, block_size: float,
       maxmiss: float, all_chromosomes: bool, output_format: str):
    """Compute f2 statistic between two populations.

    f2(A, B) measures genetic drift between populations A and B.
    """
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, [pop1, pop2], block_size, maxmiss, all_chromosomes)

        result = f2_blocks.f2(pop1, pop2)
        if result is None:
            click.echo(f"Error: Could not compute f2({pop1}, {pop2})", err=True)
            sys.exit(1)

        if output_format == "json":
            click.echo(json.dumps(_json_row(result), indent=2))
        elif output_format == "tsv":
            click.echo("pop1\tpop2\test\tse\tz\tp\tn_snps")
            click.echo(f"{pop1}\t{pop2}\t{result['est']:.6f}\t{result['se']:.6f}\t{result['z']:.3f}\t{result['p']:.6g}\t{result['n_snps']}")
        else:
            click.echo(f"f2({pop1}, {pop2}):")
            click.echo(f"  Estimate: {result['est']:.6f}")
            click.echo(f"  Std Error: {result['se']:.6f}")
            click.echo(f"  Z-score: {result['z']:.3f}  (p = {result['p']:.3g})")
            click.echo(f"  N SNPs: {result['n_snps']}")

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--pop1", required=True, callback=parse_population_list,
              help="Comma-separated populations; alone, every distinct statistic on them")
@click.option("--pop2", callback=parse_population_list, help="Comma-separated populations")
@click.option("--pop3", callback=parse_population_list, help="Comma-separated populations")
@click.option("--pop4", callback=parse_population_list, help="Comma-separated populations")
@click.option("--d-mode", is_flag=True, default=False,
              help="D-statistics instead of f4 (needs --input)")
@stat_source_options
@with_options(scan_output_options)
def qpdstat(pop1: list, pop2: list, pop3: list, pop4: list, d_mode: bool, input_prefix: str | None,
            f2_dir: str | None, block_size: float, all_chromosomes: bool, output_format: str,
            output: str | None):
    """f4- or D-statistics f(pop1, pop2; pop3, pop4) for many combinations, as admixtools2 qpdstat.

    With --pop2, --pop3 and --pop4, every combination of one population from
    each list; with --pop1 alone, every distinct statistic on those
    populations. Columns: pop1-pop4, est, se, z, p and (from genotypes) n,
    the number of SNPs used.
    """
    try:
        data = stat_source(input_prefix, f2_dir, pop1 + pop2 + pop3 + pop4)
        df = pm.qpdstat(data, pop1, pop2 or None, pop3 or None, pop4 or None, f4mode=not d_mode,
                        block_size=block_size, auto_only=not all_chromosomes,
                        pseudohaploid=_pseudohaploid())
        write_table(df, output_format, output)
    except click.ClickException:
        raise
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--pop1", required=True, callback=parse_population_list,
              help="Comma-separated targets; alone, every distinct f3 on these populations")
@click.option("--pop2", callback=parse_population_list, help="Comma-separated populations")
@click.option("--pop3", callback=parse_population_list, help="Comma-separated populations")
@click.option("--outgroup-mode", is_flag=True, default=False,
              help="With --input: do not normalize by the target's heterozygosity (outgroup f3)")
@click.option("--no-bias-correction", is_flag=True, default=False,
              help="With --input: do not subtract the target's sampling bias")
@stat_source_options
@with_options(scan_output_options)
def qp3pop(pop1: list, pop2: list, pop3: list, outgroup_mode: bool, no_bias_correction: bool,
           input_prefix: str | None, f2_dir: str | None, block_size: float, all_chromosomes: bool,
           output_format: str, output: str | None):
    """f3-statistics f3(pop1; pop2, pop3) for many combinations, as admixtools2 qp3pop.

    From genotypes, f3 is normalized by the target's heterozygosity unless
    --outgroup-mode is given; from f2 blocks it is never normalized.
    """
    if f2_dir:
        for flag, given in [("--outgroup-mode", outgroup_mode), ("--no-bias-correction", no_bias_correction)]:
            if given:
                raise click.UsageError(f"{flag} applies to genotypes (--input), not to --f2-dir")
    try:
        data = stat_source(input_prefix, f2_dir, pop1 + pop2 + pop3)
        df = pm.qp3pop(data, pop1, pop2 or None, pop3 or None, outgroupmode=outgroup_mode,
                       apply_corr=not no_bias_correction, block_size=block_size,
                       auto_only=not all_chromosomes, pseudohaploid=_pseudohaploid())
        write_table(df, output_format, output)
    except click.ClickException:
        raise
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command("outgroup-f3")
@click.option("--outgroup", required=True, help="Outgroup population (the f3 target)")
@click.option("--populations", required=True, callback=parse_population_list,
              help="Comma-separated populations: the matrix rows and columns")
@with_options(matrix_options)
@stat_source_options
@with_options(scan_output_options)
def outgroup_f3(outgroup: str, populations: list, value: str, order: str, heatmap: str | None,
                input_prefix: str | None, f2_dir: str | None, block_size: float,
                all_chromosomes: bool, output_format: str, output: str | None):
    """Matrix of outgroup f3(outgroup; A, B) for every pair of populations.

    Higher values mean more shared drift. With --heatmap, also a heatmap.
    """
    try:
        data = stat_source(input_prefix, f2_dir, [outgroup] + populations)
        m = pm.outgroup_f3_matrix(data, outgroup, populations, value=value, block_size=block_size,
                                  auto_only=not all_chromosomes, pseudohaploid=_pseudohaploid())
        write_matrix(m, order, output_format, output, heatmap, f"f3({outgroup}; X, Y) {value}",
                     None if value != "z" else 0.0)
    except click.ClickException:
        raise
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command("dstat-matrix")
@click.option("--populations", required=True, callback=parse_population_list,
              help="Comma-separated populations X, Y: the matrix rows and columns")
@click.option("--pop1", default=None, help="Fixed pop1: entries are D(pop1, pop2; X, Y)")
@click.option("--pop2", default=None, help="Fixed pop2")
@click.option("--pop3", default=None, help="Fixed pop3: entries are D(X, Y; pop3, pop4)")
@click.option("--pop4", default=None, help="Fixed pop4")
@click.option("--f4-mode", is_flag=True, default=False, help="f4 instead of D")
@with_options(matrix_options)
@stat_source_options
@with_options(scan_output_options)
def dstat_matrix(populations: list, pop1: str | None, pop2: str | None, pop3: str | None,
                 pop4: str | None, f4_mode: bool, value: str, order: str, heatmap: str | None,
                 input_prefix: str | None, f2_dir: str | None, block_size: float,
                 all_chromosomes: bool, output_format: str, output: str | None):
    """Matrix of D- (or f4-) statistics for every pair X, Y of populations.

    Give --pop3 and --pop4 for D(X, Y; pop3, pop4), or --pop1 and --pop2 for
    D(pop1, pop2; X, Y). The matrix is antisymmetric.
    """
    try:
        fixed = [p for p in (pop1, pop2, pop3, pop4) if p]
        data = stat_source(input_prefix, f2_dir, fixed + populations)
        m = pm.dstat_matrix(data, populations, pop1=pop1, pop2=pop2, pop3=pop3, pop4=pop4, value=value,
                            f4mode=f4_mode, block_size=block_size, auto_only=not all_chromosomes,
                            pseudohaploid=_pseudohaploid())
        stat = "f4" if f4_mode else "D"
        title = f"{stat}(X, Y; {pop3}, {pop4})" if pop3 else f"{stat}({pop1}, {pop2}; X, Y)"
        write_matrix(m, order, output_format, output, heatmap, f"{title} {value}",
                     0.0 if value in ("est", "z") else None)
    except click.ClickException:
        raise
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--target", required=True, help="Target population (potentially admixed)")
@click.option("--source1", required=True, help="First source population")
@click.option("--source2", required=True, help="Second source population")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def f3(target: str, source1: str, source2: str, input_prefix: str | None, f2_dir: str | None,
       block_size: float, maxmiss: float, all_chromosomes: bool, output_format: str):
    """Compute f3 statistic (admixture test).

    f3(A; B, C) tests whether A is admixed between B and C.
    A negative value indicates admixture.
    """
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, [target, source1, source2], block_size, maxmiss, all_chromosomes)

        result = f2_blocks.f3(target, source1, source2)
        if result is None:
            click.echo(f"Error: Could not compute f3({target}; {source1}, {source2})", err=True)
            sys.exit(1)

        if output_format == "json":
            click.echo(json.dumps(_json_row(result), indent=2))
        elif output_format == "tsv":
            click.echo("target\tsource1\tsource2\test\tse\tz\tp\tn_blocks")
            click.echo(f"{target}\t{source1}\t{source2}\t{result['est']:.6f}\t{result['se']:.6f}\t{result['z']:.3f}\t{result['p']:.6g}\t{result['n_blocks']}")
        else:
            click.echo(f"f3({target}; {source1}, {source2}):")
            click.echo(f"  Estimate: {result['est']:.6f}")
            click.echo(f"  Std Error: {result['se']:.6f}")
            click.echo(f"  Z-score: {result['z']:.3f}  (p = {result['p']:.3g})")
            click.echo(f"  N blocks: {result['n_blocks']}")
            if result['est'] < 0 and result['z'] < -3:
                click.echo("  => Significant admixture signal")

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--pop1", required=True, help="Population 1")
@click.option("--pop2", required=True, help="Population 2")
@click.option("--pop3", required=True, help="Population 3")
@click.option("--pop4", required=True, help="Population 4")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def f4(pop1: str, pop2: str, pop3: str, pop4: str, input_prefix: str | None, f2_dir: str | None,
       block_size: float, maxmiss: float, all_chromosomes: bool, output_format: str):
    """Compute the f4 statistic f4(pop1, pop2; pop3, pop4).

    f4(A, B; C, D) is nonzero when the A-B and C-D branches share drift,
    as gene flow between them makes them do; its sign tells which way.
    This prints f4 itself: for the D-statistic (ABBA-BABA), its
    normalised form, use qpdstat --d-mode.
    """
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, [pop1, pop2, pop3, pop4], block_size, maxmiss, all_chromosomes)

        result = f2_blocks.f4(pop1, pop2, pop3, pop4)
        if result is None:
            click.echo(f"Error: Could not compute f4({pop1}, {pop2}; {pop3}, {pop4})", err=True)
            sys.exit(1)

        if output_format == "json":
            click.echo(json.dumps(_json_row(result), indent=2))
        elif output_format == "tsv":
            click.echo("pop1\tpop2\tpop3\tpop4\test\tse\tz\tp\tn_blocks")
            click.echo(f"{pop1}\t{pop2}\t{pop3}\t{pop4}\t{result['est']:.6f}\t{result['se']:.6f}\t{result['z']:.3f}\t{result['p']:.6g}\t{result['n_blocks']}")
        else:
            click.echo(f"f4({pop1}, {pop2}; {pop3}, {pop4}):")
            click.echo(f"  Estimate: {result['est']:.6f}")
            click.echo(f"  Std Error: {result['se']:.6f}")
            click.echo(f"  Z-score: {result['z']:.3f}  (p = {result['p']:.3g})")
            click.echo(f"  N blocks: {result['n_blocks']}")
            if abs(result['z']) > 3:
                click.echo("  => Significant deviation from tree-like evolution")

    except EXPECTED_ERRORS as e:
        fail(e)
