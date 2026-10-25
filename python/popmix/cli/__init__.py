"""
Popmix CLI - Command-line interface for population genetics analysis.

Usage:
    popmix [OPTIONS] COMMAND [ARGS]...

Commands:
    info        Show dataset information
    populations List all populations with sample counts
    extract-f2  Extract F2 blocks from genotype data
    f2          Compute f2 statistic between two populations
    f3          Compute f3 statistic (admixture test)
    f4          Compute f4 statistic (D-statistic)
    qpadm       Estimate admixture weights
    qpwave      Test admixture wave count
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
    alder       ALDER weighted LD: admixture dates and test
    rolloff     ROLLOFF: admixture dates from weighted LD correlation
"""

from __future__ import annotations

import json
import math
import sys

import click

import popmix as pm


def _p_text(p: float, digits: int = 4) -> str:
    """A p-value for text output: NA where there is none (a saturated model)."""
    return "NA" if math.isnan(p) else f"{p:.{digits}f}"


def _json_number(x: float) -> float | None:
    """NaN as null, so JSON output stays valid JSON."""
    return None if math.isnan(x) else x


def parse_population_list(ctx, param, value: str) -> list:
    """Parse comma-separated population list."""
    if not value:
        return []
    return [p.strip() for p in value.split(",") if p.strip()]


def parse_snp_range(ctx, param, value: str | None) -> tuple[int, int] | None:
    """Parse SNP range in format 'start:end'."""
    if not value:
        return None
    try:
        parts = value.split(":")
        if len(parts) != 2:
            raise click.BadParameter("SNP range must be in format 'start:end'")
        return (int(parts[0]), int(parts[1]))
    except ValueError:
        raise click.BadParameter("SNP range must be integers in format 'start:end'") from None


def format_population_error(invalid_pop: str, prefix: str) -> str:
    """Format helpful error message for invalid population name."""
    lines = [f"Population '{invalid_pop}' not found in dataset."]

    try:
        suggestions = pm.suggest_populations(prefix, invalid_pop, limit=5)
        if suggestions:
            lines.append("\nDid you mean one of these?")
            for pop, count in suggestions:
                lines.append(f"  - {pop} ({count} samples)")
    except Exception:
        pass  # If suggestion fails, just skip it

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
    except Exception as e:
        click.echo(f"Warning: Could not validate populations: {e}", err=True)
        return

    missing = [p for p in populations if p not in available]
    if missing:
        # Format error message for the first missing population
        error_msg = format_population_error(missing[0], prefix)
        if len(missing) > 1:
            error_msg += f"\n\nAlso missing: {', '.join(missing[1:])}"
        raise click.ClickException(error_msg)


def _store_pseudohaploid(ctx: click.Context, param: click.Parameter, value: str) -> str:
    ctx.meta["pseudohaploid"] = value
    return value


def _pseudohaploid() -> str:
    """The --pseudohaploid rule of the running command."""
    return click.get_current_context().meta.get("pseudohaploid", "detect")


pseudohaploid_option = click.option(
    "--pseudohaploid", type=click.Choice(["detect", "admixtools2", "diploid"]), default="detect",
    expose_value=False, callback=_store_pseudohaploid,
    help="With --input: which samples count one allele per SNP. detect (default): those with no "
         "heterozygous call; admixtools2: none among the first 1,000 SNPs, which also catches some "
         "diploid genomes; diploid: none")


def f2_source_options(command):
    """Add the options that say where an f-statistic command gets its f2 blocks."""
    options = [
        click.option("--input", "-i", "input_prefix", default=None,
                     help="Input genotype file prefix; f2 blocks are extracted for the command's populations"),
        click.option("--f2-dir", default=None, type=click.Path(exists=True, file_okay=False),
                     help="Read f2 blocks saved by 'popmix extract-f2' instead of genotypes"),
        click.option("--block-size", type=float, default=0.05,
                     help="With --input: block size in Morgans (default: 0.05)"),
        click.option("--maxmiss", type=click.FloatRange(0, 1), default=0.0,
                     help="With --input: drop SNPs missing in more than this fraction of populations "
                          "(default: 0)"),
        click.option("--all-chromosomes", is_flag=True, default=False,
                     help="With --input: use all chromosomes, not only autosomes 1-22"),
        pseudohaploid_option,
    ]
    for option in reversed(options):
        command = option(command)
    return command


def load_f2_blocks(input_prefix: str | None, f2_dir: str | None, populations: list,
                   block_size: float, maxmiss: float, all_chromosomes: bool):
    """f2 blocks for POPULATIONS, as admixtools2's extract_f2() computes them.

    From genotypes, the blocks are extracted for exactly these populations,
    so with maxmiss = 0 the SNPs are those with data in all of them. To use
    the SNP set of a larger population set, as admixtools2 does with
    precomputed f2, extract it once with 'popmix extract-f2' and pass --f2-dir.
    """
    if input_prefix is None and f2_dir is None:
        raise click.UsageError("Missing --input (genotypes) or --f2-dir (saved f2 blocks)")
    if input_prefix is not None and f2_dir is not None:
        raise click.UsageError("Give --input or --f2-dir, not both")
    pops = list(dict.fromkeys(populations))
    if f2_dir is not None:
        f2_blocks = pm.read_f2(f2_dir)
        available = set(f2_blocks.populations)
        missing = [p for p in pops if p not in available]
        if missing:
            raise click.ClickException(f"Not in the f2 blocks in {f2_dir}: {', '.join(missing)}")
        return f2_blocks
    validate_populations(input_prefix, pops)
    return pm.extract_f2(input_prefix, pops, auto_only=not all_chromosomes, block_size=block_size,
                         maxmiss=maxmiss, pseudohaploid=_pseudohaploid())


@click.group()
@click.option("--threads", "-t", type=int, default=4, help="Number of threads to use (default: 4)")
@click.version_option(version=pm.__version__, prog_name="popmix")
@click.pass_context
def cli(ctx, threads: int):
    """Popmix: Population genetics F-statistics and admixture analysis toolkit."""
    ctx.ensure_object(dict)
    if threads < 1:
        raise click.BadParameter("Thread count must be at least 1")
    ctx.obj["threads"] = threads
    pm.init_threads(threads)


@cli.command()
@click.argument("prefix")
@click.option("--populations", "-p", is_flag=True, help="Show all populations with sample counts")
def info(prefix: str, populations: bool):
    """Show dataset information.

    PREFIX is the path prefix for the genotype files (without extension),
    or the full path to a VCF file (with .vcf or .vcf.gz extension).
    """
    try:
        info_dict = pm.dataset_info(prefix)
        click.echo(f"Dataset: {prefix}")
        click.echo(f"Format: {info_dict.get('format', 'unknown')}")
        click.echo(f"Samples: {info_dict.get('n_samples', 'unknown')}")
        click.echo(f"SNPs: {info_dict.get('n_snps', 'unknown')}")

        pop_list = info_dict.get("populations", [])
        click.echo(f"Populations: {len(pop_list)}")

        if populations:
            # Show full population list with counts
            pop_counts = pm.list_populations(prefix)
            for pop, count in pop_counts:
                click.echo(f"  {pop}: {count} samples")
        elif len(pop_list) <= 10:
            for pop in pop_list:
                click.echo(f"  - {pop}")
        else:
            for pop in pop_list[:5]:
                click.echo(f"  - {pop}")
            click.echo(f"  ... and {len(pop_list) - 5} more")
            click.echo("\nTip: Use --populations to show all populations with counts")
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("prefix")
@click.option("--filter", "-f", "pattern", help="Filter populations (supports wildcards like French*)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def populations(prefix: str, pattern: str | None, output_format: str):
    """List all populations in a dataset with sample counts.

    PREFIX is the path prefix for the genotype files (without extension),
    or the full path to a VCF file (with .vcf or .vcf.gz extension).

    Examples:
        popmix populations /path/to/data
        popmix populations /path/to/data.vcf.gz
        popmix populations /path/to/data --filter "French*"
        popmix populations /path/to/data --format json
    """
    import fnmatch
    import json

    try:
        pop_counts = pm.list_populations(prefix)

        # Apply filter if specified
        if pattern:
            pop_counts = [(p, c) for p, c in pop_counts if fnmatch.fnmatch(p, pattern)]

        if output_format == "json":
            data = [{"population": p, "samples": c} for p, c in pop_counts]
            click.echo(json.dumps(data, indent=2))
        elif output_format == "tsv":
            click.echo("population\tsamples")
            for pop, count in pop_counts:
                click.echo(f"{pop}\t{count}")
        else:
            click.echo(f"Populations in {prefix} ({len(pop_counts)} {'matched' if pattern else 'total'}):")
            click.echo()
            click.echo(f"{'Population':<40} Samples")
            click.echo("-" * 50)
            for pop, count in pop_counts:
                click.echo(f"{pop:<40} {count:>6}")

            if not pattern:
                click.echo()
                if len(pop_counts) > 0:
                    example_pops = [p for p, _ in pop_counts[:3]]
                    click.echo('Tip: Use --filter "pattern*" to search')
                    click.echo(f'Tip: Specify populations with: --populations "{",".join(example_pops)}"')

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command("extract-f2")
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--output", "-o", "output_dir", required=True, help="Output directory for F2 blocks")
@click.option("--snp-range", callback=parse_snp_range,
              help="Only SNPs start:end, by index (computed outside admixtools2's extract_f2 pipeline)")
@click.option("--block-size", type=float, default=0.05, help="Block size in Morgans (default: 0.05)")
@click.option("--populations", callback=parse_population_list,
              help="Comma-separated list of populations to include (default: all)")
@click.option("--maxmiss", type=click.FloatRange(0, 1), default=0.0,
              help="Drop SNPs missing in more than this fraction of populations (default: 0)")
@click.option("--auto-only/--all-chromosomes", default=True, help="Use only autosomal SNPs (default: auto-only)")
@pseudohaploid_option
def extract_f2(input_prefix: str, output_dir: str, snp_range: tuple[int, int] | None,
               block_size: float, populations: list, maxmiss: float, auto_only: bool):
    """Extract F2 blocks from genotype data, as admixtools2's extract_f2().

    Computes pairwise F2 statistics for all population pairs and saves them
    for use in downstream analyses (f3, f4, qpAdm, etc. with --f2-dir).
    """
    try:
        if not snp_range:
            if not populations:
                populations = [p for p, _ in pm.list_populations(input_prefix)]
            else:
                validate_populations(input_prefix, populations)
            click.echo(f"Streaming F2 extraction from {input_prefix}...")
            click.echo(f"Populations: {', '.join(populations)}")
            f2_blocks = pm.extract_f2(
                input_prefix,
                populations,
                auto_only=auto_only,
                block_size=block_size,
                maxmiss=maxmiss,
                pseudohaploid=_pseudohaploid(),
            )
            click.echo(f"Computed {f2_blocks.n_blocks} blocks for {f2_blocks.n_populations} populations")
        else:
            # A SNP index range needs the genotypes in memory
            click.echo(f"Reading genotypes from {input_prefix}...")
            geno = pm.read_genotypes(input_prefix, snp_range=snp_range)
            click.echo(f"Loaded {geno.n_snps} SNPs x {geno.n_samples} samples")

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

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
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
            click.echo(json.dumps(result, indent=2))
        elif output_format == "tsv":
            click.echo("pop1\tpop2\testimate\tstd_err\tz_score\tn_snps")
            click.echo(f"{pop1}\t{pop2}\t{result['estimate']:.6f}\t{result['std_err']:.6f}\t{result['z_score']:.3f}\t{int(result['n_snps'])}")
        else:
            click.echo(f"f2({pop1}, {pop2}):")
            click.echo(f"  Estimate: {result['estimate']:.6f}")
            click.echo(f"  Std Error: {result['std_err']:.6f}")
            click.echo(f"  Z-score: {result['z_score']:.3f}")
            click.echo(f"  N SNPs: {int(result['n_snps'])}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


def stat_source_options(command):
    """Options that say where a statistic scan gets its data: genotypes or saved f2 blocks."""
    options = [
        click.option("--input", "-i", "input_prefix", default=None,
                     help="Genotype prefix; statistics are computed from the genotypes, each with "
                          "the SNPs that have data in its populations (as admixtools2 qpdstat/qp3pop "
                          "on genotypes)"),
        click.option("--f2-dir", default=None, type=click.Path(exists=True, file_okay=False),
                     help="Use f2 blocks saved by 'popmix extract-f2' instead (f4 and f3 only)"),
        click.option("--block-size", type=float, default=0.05,
                     help="With --input: block size in Morgans (default: 0.05)"),
        click.option("--all-chromosomes", is_flag=True, default=False,
                     help="With --input: use all chromosomes, not only autosomes 1-22"),
        pseudohaploid_option,
    ]
    for option in reversed(options):
        command = option(command)
    return command


def stat_source(input_prefix: str | None, f2_dir: str | None, populations: list):
    """The genotype prefix, or the saved f2 blocks, for a statistic scan."""
    if (input_prefix is None) == (f2_dir is None):
        raise click.UsageError("Give one of --input (genotypes) or --f2-dir (saved f2 blocks)")
    pops = list(dict.fromkeys(populations))
    if input_prefix is not None:
        validate_populations(input_prefix, pops)
        return input_prefix
    f2_blocks = pm.read_f2(f2_dir)
    missing = [p for p in pops if p not in set(f2_blocks.populations)]
    if missing:
        raise click.ClickException(f"Not in the f2 blocks in {f2_dir}: {', '.join(missing)}")
    return f2_blocks


def write_table(df, output_format: str, output: str | None, index: bool = False) -> None:
    """Write a result table as TSV or JSON records to OUTPUT or stdout."""
    if output_format == "json":
        text = df.to_json(orient="index" if index else "records", indent=2, double_precision=15)
    else:
        text = df.to_csv(sep="\t", index=index, float_format="%.12g")
    if output:
        with open(output, "w") as f:
            f.write(text)
    else:
        click.echo(text, nl=False)


def write_matrix(matrix, order: str, output_format: str, output: str | None, heatmap: str | None,
                 title: str, center: float | None) -> None:
    """Order a population matrix, write it, and optionally draw its heatmap."""
    if order == "cluster":
        labels = pm.cluster_order(matrix)
        matrix = matrix.loc[labels, labels]
    write_table(matrix, output_format, output, index=True)
    if heatmap:
        from popmix import plot

        plot.heatmap(matrix, title=title, center=center, path=heatmap)


scan_output_options = [
    click.option("--format", "output_format", type=click.Choice(["tsv", "json"]), default="tsv",
                 help="Output format (default: tsv)"),
    click.option("--output", "-o", default=None, help="Write the table here instead of stdout"),
]


def with_options(options):
    def apply(command):
        for option in reversed(options):
            command = option(command)
        return command
    return apply


@cli.command()
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
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
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
    try:
        data = stat_source(input_prefix, f2_dir, pop1 + pop2 + pop3)
        df = pm.qp3pop(data, pop1, pop2 or None, pop3 or None, outgroupmode=outgroup_mode,
                       apply_corr=not no_bias_correction, block_size=block_size,
                       auto_only=not all_chromosomes, pseudohaploid=_pseudohaploid())
        write_table(df, output_format, output)
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


matrix_options = [
    click.option("--value", type=click.Choice(["est", "se", "z", "p"]), default="est",
                 help="Statistic in the matrix (default: est)"),
    click.option("--order", type=click.Choice(["input", "cluster"]), default="input",
                 help="Row and column order: as given, or clustered so similar populations are adjacent"),
    click.option("--heatmap", default=None, help="Also draw a heatmap to this image file (needs popmix[plot])"),
]


@cli.command("outgroup-f3")
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
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command("dstat-matrix")
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
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Genotype prefix")
@click.option("--pca-pops", callback=parse_population_list, default="",
              help="Comma-separated populations the PCA is computed on (default: all others)")
@click.option("--project-pops", callback=parse_population_list, default="",
              help="Comma-separated populations projected onto the PCA")
@click.option("--components", type=int, default=10, help="Eigenvectors reported (default: 10)")
@click.option("--outlier-iterations", type=int, default=5,
              help="Rounds of outlier removal, 0 for none (default: 5)")
@click.option("--outlier-components", type=int, default=10,
              help="Eigenvectors examined for outliers (default: 10)")
@click.option("--outlier-sigma", type=float, default=6.0,
              help="Outlier threshold in standard deviations (default: 6)")
@click.option("--lsqproject", is_flag=True, default=False,
              help="Project by least squares on observed genotypes (for samples with missing data)")
@click.option("--compat", type=click.Choice(["eigensoft-7.2.1"]), default=None,
              help="Reproduce EIGENSOFT 7.2.1 smartpca exactly (see DEVIATIONS.md)")
@click.option("--output-prefix", "-o", required=True,
              help="Writes PREFIX.coords.tsv, .eval, .tracy_widom.tsv and .outliers.tsv")
@click.option("--plot", default=None, help="Also plot PC1 against PC2 to this image file (needs popmix[plot])")
def pca(input_prefix: str, pca_pops: list, project_pops: list, components: int, outlier_iterations: int,
        outlier_components: int, outlier_sigma: float, lsqproject: bool, compat: str | None,
        output_prefix: str, plot: str | None):
    """Principal component analysis as EIGENSOFT smartpca, with projection.

    The PCA is computed on --pca-pops and the samples of --project-pops are
    projected onto it (smartpca poplistname); use --lsqproject for ancient
    samples with missing data.
    """
    try:
        validate_populations(input_prefix, pca_pops + project_pops)
        r = pm.smartpca(input_prefix, pca_pops or None, project_pops, n_components=components,
                        outlier_iterations=outlier_iterations, outlier_components=outlier_components,
                        outlier_sigma=outlier_sigma, lsqproject=lsqproject, compat=compat)
        r.coords.to_csv(f"{output_prefix}.coords.tsv", sep="\t", index=False, float_format="%.12g")
        with open(f"{output_prefix}.eval", "w") as f:
            f.write("".join(f"{x:.12g}\n" for x in r.eigenvalues))
        r.tracy_widom.to_csv(f"{output_prefix}.tracy_widom.tsv", sep="\t", index=False, float_format="%.12g")
        r.outliers.to_csv(f"{output_prefix}.outliers.tsv", sep="\t", index=False, float_format="%.12g")
        click.echo(f"{r.n_pca_samples} PCA samples, {len(r.coords) - r.n_pca_samples} projected, "
                   f"{r.n_snps_used} SNPs, {len(r.outliers)} outliers removed")
        if plot:
            from popmix import plot as plots

            plots.pca(r, path=plot)
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Genotype prefix")
@click.option("--populations", callback=parse_population_list, default="",
              help="Comma-separated populations whose samples are compared (default: all)")
@click.option("--method", type=click.Choice(["king", "read"]), default="king",
              help="king: KING-robust, for diploid calls (as PLINK 2 --make-king-table); "
                   "read: READ-style, for pseudo-haploid ancient samples (as READv2)")
@click.option("--norm", type=click.Choice(["median", "mean", "max"]), default="median",
              help="With --method read: normalize P0 by this statistic over the pairs (default: median)")
@click.option("--norm-value", type=float, default=None,
              help="With --method read: normalize P0 by this value instead")
@click.option("--window-size", type=int, default=5_000_000,
              help="With --method read: jackknife window in base pairs (default: 5000000)")
@click.option("--compat", type=click.Choice(["readv2-2.01"]), default=None,
              help="With --method read: reproduce READv2 exactly (see DEVIATIONS.md)")
@click.option("--output", "-o", default=None, help="Write the table here instead of stdout")
def kinship(input_prefix: str, populations: list, method: str, norm: str, norm_value: float | None,
            window_size: int, compat: str | None, output: str | None):
    """Kinship for every pair of samples: KING-robust or READ-style.

    KING-robust needs diploid genotype calls; for pseudo-haploid ancient
    samples use --method read, which compares them within one group and
    normalizes by the typical (by default median) pair.
    """
    try:
        validate_populations(input_prefix, populations)
        pops = populations or None
        if method == "king":
            df = pm.king(input_prefix, pops)
        else:
            df = pm.read_kinship(input_prefix, pops, norm=norm, norm_value=norm_value, window_size=window_size,
                                 compat=compat)
            click.echo(f"normalization value (P0 of unrelated pairs): {df.attrs['norm_value']:.6g}", err=True)
        write_table(df, "tsv", output)
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Genotype prefix")
@click.option("--populations", required=True, callback=parse_population_list,
              help="Comma-separated populations, one axis each")
@click.option("--projection", callback=parse_population_list, default="",
              help="Comma-separated allele counts to project each population to (default: full sizes)")
@click.option("--folded", is_flag=True, default=False, help="Fold the spectrum (minor allele counts)")
@click.option("--joint", is_flag=True, default=False,
              help="admixtools2 joint_sfs(): counts among observed alleles, no projection")
@click.option("--compat", type=click.Choice(["admixtools2-2.0.10"]), default=None,
              help="With --joint: reproduce admixtools2 exactly (see DEVIATIONS.md)")
@click.option("--all-chromosomes", is_flag=True, default=False,
              help="Use every SNP, not only autosomes (the default with --joint, as admixtools2)")
@click.option("--output", "-o", default=None, help="Write the spectrum table here instead of stdout")
def sfs(input_prefix: str, populations: list, projection: list, folded: bool, joint: bool,
        compat: str | None, all_chromosomes: bool, output: str | None):
    """Site-frequency spectrum, as dadi (projected) or admixtools2 joint_sfs().

    Writes one row per entry: the count in each population and n, the number
    of SNPs (expected, after projection). For one population, also prints S,
    pi, Watterson's theta and Tajima's D.
    """
    import itertools

    import pandas as pd

    try:
        validate_populations(input_prefix, populations)
        if joint:
            df = pm.joint_sfs(input_prefix, populations, compat=compat, auto_only=False)
        else:
            sizes = [int(x) for x in projection] or None
            spec = pm.sfs(input_prefix, populations, sizes, folded=folded, auto_only=not all_chromosomes)
            grid = list(itertools.product(*(range(d) for d in spec.data.shape)))
            df = pd.DataFrame(grid, columns=populations)
            df["n"] = spec.data.ravel()
            df["masked"] = spec.mask.ravel()
            if len(populations) == 1 and not folded:
                click.echo(f"S\t{spec.S():.6g}\npi\t{spec.pi():.6g}\n"
                           f"watterson_theta\t{spec.watterson_theta():.6g}\ntajima_d\t{spec.tajima_d():.6g}")
        if output is None:
            click.echo(df.to_csv(sep="\t", index=False, float_format="%.12g"), nl=False)
        else:
            df.to_csv(output, sep="\t", index=False, float_format="%.12g")
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Genotype prefix")
@click.option("--populations", callback=parse_population_list, default="",
              help="Comma-separated populations whose samples are used (default: all)")
@click.option("-k", "k", type=int, required=True, help="Number of ancestral populations")
@click.option("--supervised", callback=parse_population_list, default="",
              help="Comma-separated reference populations (K of them) for a supervised run")
@click.option("--fixed-p", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Projection: a .P file (SNPs x K) from an earlier run on the same SNPs")
@click.option("--cv", type=int, default=None, help="Folds of cross-validation")
@click.option("--chromosomes", callback=parse_population_list, default="",
              help="Comma-separated chromosomes to use (default: all)")
@click.option("--seed", type=int, default=1, help="Seed for the random start (default: 1)")
@click.option("--tolerance", type=float, default=1e-4,
              help="Stop when the log-likelihood improves by less than this (default: 1e-4)")
@click.option("--output-prefix", "-o", required=True,
              help="Writes PREFIX.K.Q.tsv (id, pop, Q1..QK), PREFIX.K.Q and PREFIX.K.P (ADMIXTURE's formats)")
@click.option("--plot", default=None, help="Also draw a bar plot of Q to this image file (needs popmix[plot])")
def admixture(input_prefix: str, populations: list, k: int, supervised: list, fixed_p: str | None,
              cv: int | None, chromosomes: list, seed: int, tolerance: float, output_prefix: str,
              plot: str | None):
    """Ancestry proportions with K ancestral populations, as ADMIXTURE.

    Unsupervised by default; --supervised fixes the reference populations'
    samples, --fixed-p projects new samples onto an earlier run's P.
    """
    import numpy as np

    try:
        validate_populations(input_prefix, populations + supervised)
        p = None if fixed_p is None else np.loadtxt(fixed_p, ndmin=2)
        r = pm.admixture(input_prefix, k, populations or None, supervised=supervised or None, fixed_p=p, cv=cv,
                         chromosomes=chromosomes or None, seed=seed, tolerance=tolerance)
        r.Q.to_csv(f"{output_prefix}.{k}.Q.tsv", sep="\t", index=False, float_format="%.12g")
        np.savetxt(f"{output_prefix}.{k}.Q", r.q_matrix(), fmt="%.6f")
        np.savetxt(f"{output_prefix}.{k}.P", r.P, fmt="%.6f")
        click.echo(f"Loglikelihood: {r.loglik:.6f} ({r.iterations} iterations"
                   f"{'' if r.converged else ', not converged'})")
        if r.cv_error is not None:
            click.echo(f"CV error (K={k}): {r.cv_error:.5f}")
        if plot:
            from popmix import plot as plots

            plots.admixture(r, path=plot)
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Genotype prefix (with a genetic map)")
@click.option("--test", "test_pop", required=True, help="Test (possibly admixed) population")
@click.option("--refs", required=True, callback=parse_population_list,
              help="One or two comma-separated reference populations")
@click.option("--mindis", type=float, default=None,
              help="Fit from this distance (Morgans) instead of where LD stops being correlated")
@click.option("--binsize", type=float, default=0.0005, help="Distance bin in Morgans (default: 0.0005)")
@click.option("--maxdis", type=float, default=0.5, help="Largest distance fit, Morgans (default: 0.5)")
@click.option("--mincount", type=int, default=4, help="Test samples called at a SNP (default: 4)")
@click.option("--compat", type=click.Choice(["alder-1.03"]), default=None,
              help="Reproduce ALDER 1.03 exactly (see DEVIATIONS.md)")
@click.option("--output-prefix", "-o", default=None,
              help="Also write PREFIX.curve.tsv and PREFIX.fits.tsv")
@click.option("--plot", default=None, help="Also plot the curve and fit to this image file (needs popmix[plot])")
def alder(input_prefix: str, test_pop: str, refs: list, mindis: float | None, binsize: float, maxdis: float,
          mincount: int, compat: str | None, output_prefix: str | None, plot: str | None):
    """ALDER weighted LD: date admixture in TEST and, with two references, test for it."""
    try:
        validate_populations(input_prefix, [test_pop] + refs)
        r = pm.alder(input_prefix, test_pop, refs, mindis=mindis, binsize=binsize, maxdis=maxdis,
                     mincount=mincount, compat=compat)
        for ref, stop in r.ld_corr_stops.items():
            if stop is not None:
                click.echo(f"LD correlated with {ref} to {100 * stop:.2f} cM")
        f = r.curve.fit
        click.echo(f"fit from {f['mindis_cm']:.2f} cM: decay {f['decay']:.2f} +/- {f['decay_se']:.2f} "
                   f"generations, amplitude {f['amp_exp']:.6g} +/- {f['amp_exp_se']:.6g}")
        if r.test is not None:
            t = r.test
            verdict = "SUCCEEDS" if t["succeeds"] else "FAILS"
            click.echo(f"admixture test {verdict} (z = {t['z']:.2f}, p = {t['p']:.2g})"
                       f"{'' if t['consistent'] else '; warning: decay rates inconsistent'}")
        if r.mixture_fraction is not None:
            m, se = r.mixture_fraction
            click.echo(f"mixture fraction lower bound: {100 * m:.1f}% +/- {100 * se:.1f}%")
        if output_prefix:
            r.curve.curve.to_csv(f"{output_prefix}.curve.tsv", sep="\t", index=False, float_format="%.12g")
            r.curve.fits.to_csv(f"{output_prefix}.fits.tsv", sep="\t", index=False, float_format="%.12g")
        if plot:
            from popmix import plot as plots

            plots.alder_curve(r, path=plot)
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Genotype prefix (with a genetic map)")
@click.option("--admixed", required=True, help="Admixed population")
@click.option("--refs", required=True, callback=parse_population_list, help="The two comma-separated references")
@click.option("--binsize", type=float, default=0.0005, help="Distance bin in Morgans (default: 0.0005)")
@click.option("--maxdis", type=float, default=0.05, help="Largest distance in Morgans (default: 0.05)")
@click.option("--runmode", type=click.Choice(["0", "1"]), default="1", help="Weights: 1 pA - pB (default), 0 scaled")
@click.option("--chithresh", type=float, default=6.0,
              help="Leave out pairs in LD in a reference, chi-square at least this (default: 6)")
@click.option("--chromosomes", callback=parse_population_list, default="",
              help="Comma-separated chromosomes to use (default: all autosomes)")
@click.option("--output-prefix", "-o", default=None, help="Also write PREFIX.bins.tsv")
def rolloff(input_prefix: str, admixed: str, refs: list, binsize: float, maxdis: float, runmode: str,
            chithresh: float, chromosomes: list, output_prefix: str | None):
    """ROLLOFF: date admixture from the decay of weighted LD correlation (ADMIXTOOLS)."""
    try:
        validate_populations(input_prefix, [admixed] + refs)
        r = pm.rolloff(input_prefix, admixed, refs, binsize=binsize, maxdis=maxdis, runmode=int(runmode),
                       chithresh=chithresh, chromosomes=chromosomes or None)
        m, se = r.decay_jackknife
        click.echo(f"decay {r.decay:.3f} generations; jackknife {m:.3f} +/- {se:.3f} "
                   f"({len(r.chromosomes)} chromosomes, {r.n_snps_used} SNPs)")
        if output_prefix:
            r.bins.to_csv(f"{output_prefix}.bins.tsv", sep="\t", index=False, float_format="%.12g")
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
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
            click.echo(json.dumps(result, indent=2))
        elif output_format == "tsv":
            click.echo("target\tsource1\tsource2\testimate\tstd_err\tz_score\tn_blocks")
            click.echo(f"{target}\t{source1}\t{source2}\t{result['estimate']:.6f}\t{result['std_err']:.6f}\t{result['z_score']:.3f}\t{int(result['n_blocks'])}")
        else:
            click.echo(f"f3({target}; {source1}, {source2}):")
            click.echo(f"  Estimate: {result['estimate']:.6f}")
            click.echo(f"  Std Error: {result['std_err']:.6f}")
            click.echo(f"  Z-score: {result['z_score']:.3f}")
            click.echo(f"  N blocks: {int(result['n_blocks'])}")
            if result['estimate'] < 0 and result['z_score'] < -3:
                click.echo("  => Significant admixture signal")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--pop1", required=True, help="Population 1")
@click.option("--pop2", required=True, help="Population 2")
@click.option("--pop3", required=True, help="Population 3")
@click.option("--pop4", required=True, help="Population 4")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def f4(pop1: str, pop2: str, pop3: str, pop4: str, input_prefix: str | None, f2_dir: str | None,
       block_size: float, maxmiss: float, all_chromosomes: bool, output_format: str):
    """Compute f4 statistic (D-statistic / ABBA-BABA test).

    f4(A, B; C, D) tests for gene flow between (A,B) and (C,D) branches.
    """
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, [pop1, pop2, pop3, pop4], block_size, maxmiss, all_chromosomes)

        result = f2_blocks.f4(pop1, pop2, pop3, pop4)
        if result is None:
            click.echo(f"Error: Could not compute f4({pop1}, {pop2}; {pop3}, {pop4})", err=True)
            sys.exit(1)

        if output_format == "json":
            click.echo(json.dumps(result, indent=2))
        elif output_format == "tsv":
            click.echo("pop1\tpop2\tpop3\tpop4\testimate\tstd_err\tz_score\tn_blocks")
            click.echo(f"{pop1}\t{pop2}\t{pop3}\t{pop4}\t{result['estimate']:.6f}\t{result['std_err']:.6f}\t{result['z_score']:.3f}\t{int(result['n_blocks'])}")
        else:
            click.echo(f"f4({pop1}, {pop2}; {pop3}, {pop4}):")
            click.echo(f"  Estimate: {result['estimate']:.6f}")
            click.echo(f"  Std Error: {result['std_err']:.6f}")
            click.echo(f"  Z-score: {result['z_score']:.3f}")
            click.echo(f"  N blocks: {int(result['n_blocks'])}")
            if abs(result['z_score']) > 3:
                click.echo("  => Significant deviation from tree-like evolution")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
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

        result = pm.qpadm(f2_blocks, target, left, right, constrained=constrained)

        if output_format == "json":
            output = {
                "target": result.target,
                "sources": result.sources,
                "weights": dict(zip(result.sources, result.weights)),
                "std_errors": dict(zip(result.sources, result.weight_std_errors)),
                "chi_sq": result.chi_sq,
                "df": result.df,
                "p_value": _json_number(result.p_value),
                "feasible": result.is_feasible(),
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("target\tsource\tweight\tchi_sq\tdf\tp_value")
            for source, weight in zip(result.sources, result.weights):
                click.echo(f"{target}\t{source}\t{weight:.4f}\t{result.chi_sq:.4f}\t{result.df}\t{_p_text(result.p_value)}")
        else:
            click.echo(f"qpAdm: {target} = mix of {', '.join(left)}")
            click.echo(f"Outgroups: {', '.join(right)}")
            click.echo()
            click.echo("Admixture weights:")
            for source, weight in zip(result.sources, result.weights):
                click.echo(f"  {source}: {weight:.4f} ({weight*100:.1f}%)")
            click.echo()
            click.echo(f"Model fit: chi-sq={result.chi_sq:.4f}, df={result.df}, p={_p_text(result.p_value)}")
            click.echo(f"Feasible: {result.is_feasible()}")
            if math.isnan(result.p_value):
                click.echo("=> Model not tested: it has 0 degrees of freedom (as many right populations as"
                           " left), so the data cannot reject it")
            elif result.is_good_fit():
                click.echo("=> Good model fit (p > 0.05)")
            else:
                click.echo("=> Poor model fit (p <= 0.05)")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
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
                "left": result.left_pops,
                "right": result.right_pops,
                "estimated_rank": result.estimated_rank,
                "rank_tests": [
                    {"rank": t.rank, "chi_sq": t.chi_sq, "df": t.df, "p_value": _json_number(t.p_value)}
                    for t in result.rank_tests
                ],
                "is_clade": result.is_clade(),
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("rank\tchi_sq\tdf\tp_value")
            for t in result.rank_tests:
                click.echo(f"{t.rank}\t{t.chi_sq:.4f}\t{t.df}\t{_p_text(t.p_value, 6)}")
        else:
            click.echo(f"qpWave: {', '.join(left)}")
            click.echo(f"Outgroups: {', '.join(right)}")
            click.echo()
            click.echo("Rank tests:")
            click.echo("  Rank  Chi-sq     df    P-value")
            for t in result.rank_tests:
                sig = "*" if t.p_value < 0.05 else " "
                click.echo(f"  {t.rank:4d}  {t.chi_sq:8.4f}  {t.df:4d}  {_p_text(t.p_value, 6):>8} {sig}")
            click.echo()
            click.echo(f"Estimated rank: {result.estimated_rank}")
            if result.is_clade():
                click.echo("=> Populations form a clade (rank 0 not rejected)")
            else:
                click.echo(f"=> {result.estimated_rank + 1} independent ancestry streams")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix (or full path for VCF)")
@click.option("--output", "-o", "output_prefix", required=True, help="Output genotype file prefix")
@click.option("--input-format", "-if", "input_format", default=None,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap", "vcf"]),
              help="Input format (auto-detected if not specified)")
@click.option("--output-format", "-of", "output_format", required=True,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap"]),
              help="Output format")
def convertf(input_prefix: str, output_prefix: str, input_format: str | None, output_format: str):
    """Convert between genotype file formats.

    Supports conversion between 6 input formats:
    - vcf: VCF format (.vcf, .vcf.gz) - pass full path
    - ancestrymap: Sparse text (.ancestrymapgeno, .snp, .ind)
    - eigenstrat: Dense text (.geno, .snp, .ind)
    - packedancestrymap: Binary EIGENSTRAT (.geno, .snp, .ind)
    - packedped/plink: Binary PLINK (.bed, .bim, .fam)
    - ped: Text PLINK (.ped, .map, .pedind)

    Note: VCF output is not supported - use plink, eigenstrat, or packedancestrymap.
    """
    try:
        # Detect input format if not specified
        if input_format is None:
            detected = pm.detect_format(input_prefix)
            click.echo(f"Auto-detected input format: {detected}")
        else:
            click.echo(f"Input format: {input_format}")

        click.echo(f"Output format: {output_format}")
        click.echo(f"Converting {input_prefix} -> {output_prefix}...")

        result = pm.convertf(
            input_prefix=input_prefix,
            output_prefix=output_prefix,
            output_format=output_format,
            input_format=input_format,
        )

        click.echo(f"Converted {result['n_snps']} SNPs x {result['n_samples']} samples")
        click.echo(f"Output files written to {output_prefix}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--file-a", "-a", required=True, help="First input genotype file (prefix or full VCF path)")
@click.option("--file-b", "-b", required=True, help="Second input genotype file (prefix or full VCF path)")
@click.option("--output", "-o", required=True, help="Output genotype file prefix")
@click.option("--output-format", "-of", required=True,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap"]),
              help="Output format")
@click.option("--input-format-a", "-ifa", default=None,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap", "vcf"]),
              help="Format for file A (auto-detected if not specified)")
@click.option("--input-format-b", "-ifb", default=None,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap", "vcf"]),
              help="Format for file B (auto-detected if not specified)")
def merge(file_a: str, file_b: str, output: str, output_format: str,
          input_format_a: str | None, input_format_b: str | None):
    """Merge two genotype datasets.

    Combines samples from two genotype files, keeping only SNPs that
    exist in both datasets, paired by position (and by name where several
    share one). Alleles are reconciled as AdmixTools mergeit does: genotypes
    are recoded when the alleles are swapped, a SNP on the other strand is
    complemented, the unknown allele of a monomorphic SNP is filled in, and
    A/T, C/G and third-allele SNPs are left out and counted.

    Supports VCF input (pass full path with .vcf or .vcf.gz extension).
    VCF output is not supported.
    """
    try:
        # Detect formats if not specified
        if input_format_a is None:
            detected_a = pm.detect_format(file_a)
            click.echo(f"File A format (auto-detected): {detected_a}")
        else:
            click.echo(f"File A format: {input_format_a}")

        if input_format_b is None:
            detected_b = pm.detect_format(file_b)
            click.echo(f"File B format (auto-detected): {detected_b}")
        else:
            click.echo(f"File B format: {input_format_b}")

        click.echo(f"Output format: {output_format}")
        click.echo(f"Merging {file_a} + {file_b} -> {output}...")

        result = pm.merge(
            file_a=file_a,
            file_b=file_b,
            output=output,
            output_format=output_format,
            input_format_a=input_format_a,
            input_format_b=input_format_b,
        )

        click.echo()
        click.echo("Merge complete:")
        click.echo(f"  SNPs in file A: {result['snps_a']}")
        click.echo(f"  SNPs in file B: {result['snps_b']}")
        click.echo(f"  SNPs in common: {result['snps_merged']}")
        click.echo(f"  Samples in file A: {result['samples_a']}")
        click.echo(f"  Samples in file B: {result['samples_b']}")
        click.echo(f"  Total samples: {result['samples_merged']}")

        if result['flipped_snps'] > 0:
            click.echo(f"  Alleles swapped (recoded): {result['flipped_snps']}")
        if result['other_strand_snps'] > 0:
            click.echo(f"  On the other strand (complemented): {result['other_strand_snps']}")
        dropped = [
            (result['strand_ambiguous_snps'], "A/T or C/G (strand unknown)"),
            (result['snps_mismatch'], "allele mismatch"),
            (result['duplicate_position_snps'], "several SNPs of file B at the position, none with the name"),
        ]
        for count, reason in dropped:
            if count > 0:
                click.echo(f"  Left out, {reason}: {count}")

        click.echo(f"\nOutput files written to {output}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--test-pop", required=True, help="Test (admixed) population")
@click.option("--source1", required=True, help="First source population")
@click.option("--source2", required=True, help="Second source population")
@click.option("--bin-width", type=float, default=0.1, help="Distance bin width in cM (default: 0.1)")
@click.option("--max-distance", type=float, default=100.0, help="Maximum genetic distance in cM (default: 100)")
@click.option("--min-distance", type=float, default=0.45,
              help="Fit bins from this distance in cM, DATES lovalfit (default: 0.45)")
@click.option("--qbin", type=click.IntRange(min=0), default=0,
              help="DATES qbin: 0 visits every SNP pair; N > 0 uses an FFT mesh of N cells "
                   "per bin, faster and approximate (default: 0)")
@click.option("--compat", type=click.Choice(["dates-v753"]), default=None,
              help="Reproduce DATES v753's known quirks instead of correcting them, "
                   "for output identical to DATES (see DEVIATIONS.md)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def dates(input_prefix: str, test_pop: str, source1: str, source2: str,
          bin_width: float, max_distance: float, min_distance: float, qbin: int,
          compat: str | None, output_format: str):
    """Estimate admixture date using DATES algorithm.

    DATES estimates the time since admixture by fitting an exponential decay
    to weighted linkage disequilibrium as a function of genetic distance.
    """
    # Validate populations early
    validate_populations(input_prefix, [test_pop, source1, source2])

    try:
        click.echo(f"Loading genotypes from {input_prefix}...", err=True)
        geno = pm.read_genotypes(input_prefix)
        click.echo(f"Loaded {geno.n_snps} SNPs x {geno.n_samples} samples", err=True)

        click.echo(f"Running DATES: {test_pop} = mix of {source1} and {source2}...", err=True)
        result = pm.dates(
            geno,
            test_pop=test_pop,
            source1=source1,
            source2=source2,
            bin_width_cm=bin_width,
            max_distance_cm=max_distance,
            min_distance_cm=min_distance,
            qbin=qbin,
            compat=compat,
        )

        if output_format == "json":
            output = {
                "test_pop": test_pop,
                "source1": source1,
                "source2": source2,
                "amplitude": result.fit.amplitude,
                "generations": result.generations,
                "generations_se": result.generations_se,
                "years_bp": result.years_bp,
                "years_bp_se": result.years_bp_se,
                "n_snps": result.n_snps,
                "n_bins": len(result.ld_bins),
                "qbin": qbin,
                "compat": compat,
                "fit": {
                    "amplitude": result.fit.amplitude,
                    "decay_rate": result.fit.decay_rate,
                    "r_squared": result.fit.r_squared,
                },
                "ld_bins": [
                    {
                        "distance_cm": b.distance_cm,
                        "n_pairs": b.n_pairs,
                        "correlation": b.correlation,
                        "regression_coef": b.regression_coef,
                        "weighted_ld": b.weighted_ld,
                    }
                    for b in result.ld_bins[:20]  # First 20 bins
                ],
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("distance_cm\tn_pairs\tcorrelation\tweighted_ld")
            for b in result.ld_bins:
                click.echo(f"{b.distance_cm:.1f}\t{b.n_pairs}\t{b.correlation:.6f}\t{b.weighted_ld:.6e}")
        else:
            click.echo(f"\nDATES Results: {test_pop} = mix of {source1} and {source2}")
            click.echo("=" * 60)
            click.echo(f"Amplitude (fit): {result.fit.amplitude:.6f}")
            click.echo(f"Estimated generations: {result.generations:.1f} ± {result.generations_se:.1f}")
            click.echo(f"Estimated years BP: {result.years_bp:.0f} ± {result.years_bp_se:.0f} (29 years/gen)")
            click.echo(f"Fit R²: {result.fit.r_squared:.4f}")
            click.echo(f"SNPs used: {result.n_snps}")
            click.echo(f"LD bins: {len(result.ld_bins)}")
            click.echo()
            click.echo("LD Decay (first 10 bins):")
            click.echo("  Dist(cM)  Pairs      Correlation")
            for b in result.ld_bins[:10]:
                click.echo(f"  {b.distance_cm:6.1f}    {b.n_pairs:9d}  {b.correlation:+.6f}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


def read_graph_file(path: str) -> list:
    """Read an admixture graph: one edge per line, "from to" (whitespace-
    separated). Blank lines and text after "#" are ignored."""
    edges = []
    with open(path) as f:
        for n, line in enumerate(f, 1):
            fields = line.split("#", 1)[0].split()
            if not fields:
                continue
            if len(fields) != 2:
                raise click.BadParameter(f"{path}:{n}: expected 'from to', got {line.strip()!r}",
                                         param_hint="--graph")
            edges.append((fields[0], fields[1]))
    if not edges:
        raise click.BadParameter(f"{path} has no edges", param_hint="--graph")
    return edges


def graph_leaves(edges: list) -> list:
    """Leaves (nodes without outgoing edges) in order of first appearance."""
    order, sources = [], set()
    for a, b in edges:
        sources.add(a)
        for v in (a, b):
            if v not in order:
                order.append(v)
    return [v for v in order if v not in sources]


@cli.command()
@click.option("--graph", "graph_file", required=True, type=click.Path(exists=True, dir_okay=False),
              help='Graph file: one edge per line, "from to"; the leaves are the populations')
@click.option("--admix-weights", default=None,
              help="Fit at these admixture weights (comma-separated, one per admixture node) "
                   "instead of optimizing them")
@click.option("--numstart", type=click.IntRange(min=1), default=10,
              help="Random starts for the admixture weights (default: 10)")
@click.option("--seed", type=int, default=1, help="Seed for the random starts (default: 1)")
@click.option("--diag", type=float, default=1e-4, help="Edge-length regularization (default: 1e-4)")
@click.option("--diag-f3", type=float, default=1e-5, help="f3 covariance regularization (default: 1e-5)")
@click.option("--unconstrained", is_flag=True, default=False,
              help="Allow negative drift lengths and weights outside [0, 1]")
@click.option("--plot", default=None, help="Also draw the fitted graph to this image file (needs popmix[plot])")
@click.option("--plot-fit", default=None,
              help="Also plot each f3's residual z to this image file (needs popmix[plot])")
@click.option("--dot", default=None, help="Also write the fitted graph in Graphviz DOT format to this file")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpgraph(graph_file: str, admix_weights: str | None, numstart: int, seed: int, diag: float,
            diag_f3: float, unconstrained: bool, plot: str | None, plot_fit: str | None, dot: str | None,
            input_prefix: str | None, f2_dir: str | None, block_size: float, maxmiss: float,
            all_chromosomes: bool, output_format: str):
    """Fit an admixture graph, as admixtools2 qpgraph does.

    Fits drift lengths and admixture weights to f3(base; i, j) for all pairs
    of the graph's leaves; the first leaf in the file is the base. A node
    with two incoming edges is an admixture node, and --admix-weights gives
    the weight of its incoming edge whose source node appears first.
    """
    edges = read_graph_file(graph_file)
    pops = graph_leaves(edges)
    weights = None
    if admix_weights is not None:
        try:
            weights = [float(w) for w in admix_weights.split(",")]
        except ValueError:
            raise click.BadParameter(f"not a list of numbers: {admix_weights!r}",
                                     param_hint="--admix-weights") from None

    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, pops, block_size, maxmiss, all_chromosomes)
        result = pm.qpgraph(f2_blocks, edges, admix_weights=weights, numstart=numstart, seed=seed,
                            diag=diag, diag_f3=diag_f3, constrained=not unconstrained)
        if dot:
            pm.write_dot(result, path=dot)
        if plot or plot_fit:
            from popmix import plot as plots

            if plot:
                plots.graph(result, path=plot)
            if plot_fit:
                plots.qpgraph_fit(result, path=plot_fit)

        if output_format == "json":
            output = {
                "pops": result.pops,
                "score": result.score,
                "edges": [{"from": e.from_node, "to": e.to_node, "type": e.kind, "weight": e.weight,
                           "low": e.low, "high": e.high} for e in result.edges],
                "f3": [{"pop1": f.pop1, "pop2": f.pop2, "pop3": f.pop3, "est": f.est, "se": f.se,
                        "fit": f.fit, "diff": f.diff, "z": f.z, "p": f.p} for f in result.f3],
                "opt": [{"start": o.start, "weights": o.weights, "score": o.score,
                         "evaluations": o.evaluations, "convergence": o.convergence} for o in result.opt],
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("from\tto\ttype\tweight\tlow\thigh")
            for e in result.edges:
                low = "" if e.low is None else f"{e.low:.6g}"
                high = "" if e.high is None else f"{e.high:.6g}"
                click.echo(f"{e.from_node}\t{e.to_node}\t{e.kind}\t{e.weight:.6g}\t{low}\t{high}")
        else:
            worst = max(result.f3, key=lambda f: abs(f.z))
            click.echo(f"qpGraph: {len(pops)} populations, {len(edges)} edges")
            click.echo(f"f3 base population: {result.pops[0]}")
            click.echo(f"Score: {result.score:.4f}")
            click.echo(f"Worst f3 residual: f3({worst.pop1}; {worst.pop2}, {worst.pop3}) "
                       f"z = {worst.z:.2f}")
            click.echo()
            click.echo(f"  {'from':<16} {'to':<16} {'type':<6} weight")
            for e in result.edges:
                click.echo(f"  {e.from_node:<16} {e.to_node:<16} {e.kind:<6} {e.weight:.6f}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command("find-graphs")
@click.option("--populations", required=True, callback=parse_population_list,
              help="The graph's leaves (comma-separated); the first is the f3 base of every fit")
@click.option("--nadmix", type=click.IntRange(min=0), default=0, help="Number of admixture nodes (default: 0)")
@click.option("--outgroup", default=None, help="Only graphs whose root has this population as a direct child")
@click.option("--method", type=click.Choice(["stepwise", "exhaustive", "evolutionary"]), default="stepwise",
              help="stepwise (default): add populations, then admixture edges, one at a time, keeping the "
                   "--beam-width best graphs; exhaustive: every graph, as admixtools2 generate_all_graphs(); "
                   "evolutionary: admixtools2 find_graphs()'s search")
@click.option("--stop-gen", type=click.IntRange(min=1), default=100,
              help="evolutionary: generations at most (default: 100)")
@click.option("--stop-gen2", type=click.IntRange(min=1), default=15,
              help="evolutionary: stop after this many generations without improvement (default: 15)")
@click.option("--max-admix", type=click.IntRange(min=0), default=None,
              help="evolutionary: allow up to this many admixture nodes (default: --nadmix)")
@click.option("--plot", default=None, help="Also draw the best graph to this image file (needs popmix[plot])")
@click.option("--beam-width", type=click.IntRange(min=1), default=20,
              help="Graphs kept after each stepwise step (default: 20)")
@click.option("--top", type=click.IntRange(min=1), default=5, help="Graphs to report (default: 5)")
@click.option("--output-graph", type=click.Path(dir_okay=False, writable=True), default=None,
              help="Write the best graph here, in the format 'popmix qpgraph --graph' reads")
@click.option("--numstart", type=click.IntRange(min=1), default=None,
              help="Random starts for each graph's admixture weights (default: 10; evolutionary: 1, as admixtools2)")
@click.option("--seed", type=int, default=1, help="Seed for the random starts (default: 1)")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json"]), default="text",
              help="Output format")
def find_graphs(populations: list, nadmix: int, outgroup: str | None, method: str, stop_gen: int, stop_gen2: int,
                max_admix: int | None, plot: str | None, beam_width: int, top: int, output_graph: str | None,
                numstart: int | None, seed: int, input_prefix: str | None, f2_dir: str | None,
                block_size: float, maxmiss: float, all_chromosomes: bool, output_format: str):
    """Search for the admixture graphs that best fit the data.

    Each candidate graph has POPULATIONS as leaves and NADMIX admixture nodes
    and is fitted as 'popmix qpgraph' fits it.
    """
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, populations, block_size, maxmiss, all_chromosomes)
        if method == "evolutionary":
            # admixtools2 find_graphs fits each candidate from one random start
            graphs, n_evaluated = pm.find_graphs(f2_blocks, populations, nadmix=nadmix, outgroup=outgroup,
                                                 method=method, numstart=numstart or 1, seed=seed,
                                                 stop_gen=stop_gen, stop_gen2=stop_gen2, max_admix=max_admix)
        else:
            graphs, n_evaluated = pm.find_graphs(f2_blocks, populations, nadmix=nadmix, outgroup=outgroup,
                                                 method=method, beam_width=beam_width,
                                                 numstart=numstart or 10, seed=seed)
        if plot:
            from popmix import plot as plots

            plots.graph(graphs[0].edges, path=plot, title=f"score {graphs[0].score:.4f}")
        if output_graph is not None:
            with open(output_graph, "w") as f:
                f.write(f"# popmix find-graphs: score {graphs[0].score!r}\n")
                for a, b in graphs[0].edges:
                    f.write(f"{a}\t{b}\n")

        shown = graphs[:top]
        if output_format == "json":
            output = {
                "n_evaluated": n_evaluated,
                "graphs": [{"score": g.score, "edges": g.edges} for g in shown],
            }
            click.echo(json.dumps(output, indent=2))
        else:
            click.echo(f"find-graphs: {len(populations)} populations, {nadmix} admixture, {method} search")
            click.echo(f"Graphs fitted: {n_evaluated}")
            for rank, g in enumerate(shown, 1):
                click.echo()
                click.echo(f"#{rank}  score {g.score:.4f}")
                for a, b in g.edges:
                    click.echo(f"  {a}\t{b}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


def main():
    """Entry point for the popmix CLI."""
    cli()


if __name__ == "__main__":
    main()
