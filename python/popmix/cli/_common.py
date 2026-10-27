"""Options, sources and output shared by the commands."""

from __future__ import annotations

import math

import click

import popmix as pm

from ._population_checks import validate_populations


def _p_text(p: float, digits: int = 4) -> str:
    """A p-value for text output: NA where there is none (a saturated model)."""
    return "NA" if math.isnan(p) else f"{p:.{digits}f}"


def _json_number(x: float) -> float | None:
    """NaN as null, so JSON output stays valid JSON."""
    return None if math.isnan(x) else x


def _json_row(row: dict) -> dict:
    """A statistic's dict for JSON: NaN (z without an SE) as null."""
    return {k: _json_number(v) if isinstance(v, float) else v for k, v in row.items()}


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


def write_columns(columns: dict, output: str | None) -> None:
    """Write equal-length columns as TSV, exactly as write_table writes a DataFrame of them.

    For tables built straight from the Rust results, without importing pandas
    (a third of a fast command's time).
    """
    import csv
    import io

    import numpy as np

    def cell(x) -> str:
        if isinstance(x, (float, np.floating)):
            return "" if math.isnan(x) else f"{x:.12g}"
        return str(x)

    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(list(columns))
    for row in zip(*columns.values()):
        writer.writerow([cell(x) for x in row])
    text = buf.getvalue()
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


matrix_options = [
    click.option("--value", type=click.Choice(["est", "se", "z", "p"]), default="est",
                 help="Statistic in the matrix (default: est)"),
    click.option("--order", type=click.Choice(["input", "cluster"]), default="input",
                 help="Row and column order: as given, or clustered so similar populations are adjacent"),
    click.option("--heatmap", default=None, help="Also draw a heatmap to this image file (needs popmix[plot])"),
]
