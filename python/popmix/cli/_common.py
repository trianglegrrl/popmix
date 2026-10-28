"""Options, sources and output shared by the commands."""

from __future__ import annotations

import math

import click

import popmix as pm
from popmix import _core

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
    ctx.meta.setdefault("popmix.params", {})["pseudohaploid"] = value
    return value


def _lines(path: str) -> list[list[str]]:
    """The whitespace-separated fields of each line of PATH; blank lines and
    # comments are skipped."""
    with open(path) as f:
        return [fields for fields in (line.split() for line in f) if fields and not fields[0].startswith("#")]


def _snp_ids(path: str) -> list[str]:
    """SNP IDs, the first field of each line of PATH (a .bim file's second,
    so .snp and .bim files work too)."""
    return [fields[1] if path.endswith(".bim") else fields[0] for fields in _lines(path)]


def _sample_ids(path: str) -> list[str]:
    """Sample IDs, the first field of each line of PATH (a .fam file's
    second, so .ind and .fam files work too)."""
    return [fields[1] if path.endswith(".fam") else fields[0] for fields in _lines(path)]


def _pairs(path: str) -> dict[str, str]:
    """NAME NEW per line of PATH."""
    out = {}
    for fields in _lines(path):
        if len(fields) < 2:
            raise click.BadParameter(f"{path}: each line needs a name and its new population, got {fields[0]!r}")
        out[fields[0]] = fields[1]
    return out


def _store_sample_selection(ctx: click.Context, param: click.Parameter, value):
    ctx.meta.setdefault("popmix.params", {})[param.name] = value
    if value is not None:
        read = _pairs if param.name.startswith("relabel") else _sample_ids
        ctx.meta.setdefault("sample_selection", {})[param.name] = read(value)
    return value


def sample_selection() -> dict:
    """The sample selection of the running command, as ``select_samples`` keywords."""
    return dict(click.get_current_context().meta.get("sample_selection", {}))


def _sample_selection_options(exposed: bool) -> list:
    path = click.Path(exists=True, dir_okay=False)
    kw = {"default": None, "type": path} if exposed else \
        {"default": None, "type": path, "expose_value": False, "callback": _store_sample_selection}
    return [
        click.option("--keep-inds", **kw,
                     help="Use only the sample IDs in this file (one per line, or a .ind/.fam file)"),
        click.option("--exclude-inds", **kw,
                     help="Never use the sample IDs in this file (one per line, or a .ind/.fam file)"),
        click.option("--relabel", **kw,
                     help="Give samples new populations: a file of ID POPULATION lines"),
        click.option("--relabel-pops", **kw,
                     help="Rename populations: a file of OLD NEW lines (after --relabel)"),
    ]


sample_selection_options = _sample_selection_options(exposed=False)


def _store_snp_filter(ctx: click.Context, param: click.Parameter, value):
    ctx.meta.setdefault("popmix.params", {})[param.name] = value
    filters = ctx.meta.setdefault("snp_filters", {})
    if param.name == "transversions_only" and value:
        filters["transitions"] = False
    elif param.name == "transitions_only" and value:
        filters["transversions"] = False
    elif param.name in ("minmaf", "maxmaf") and value != param.default:
        filters[param.name] = value
    elif param.name in ("keep_snps", "exclude_snps") and value is not None:
        filters[param.name.replace("_", "")] = _snp_ids(value)
    if filters.get("transitions") is False and filters.get("transversions") is False:
        raise click.UsageError("--transversions-only and --transitions-only leave no SNPs")
    return value


def snp_filters() -> dict:
    """The SNP filters of the running command, as ``extract_f2`` keywords."""
    return dict(click.get_current_context().meta.get("snp_filters", {}))


snp_filter_options = [
    click.option("--transversions-only", is_flag=True, default=False, expose_value=False,
                 callback=_store_snp_filter,
                 help="Use only transversions (drop A<->G and C<->T), for ancient DNA without UDG treatment"),
    click.option("--transitions-only", is_flag=True, default=False, expose_value=False, callback=_store_snp_filter,
                 help="Use only transitions"),
    click.option("--minmaf", type=click.FloatRange(0, 0.5), default=0.0, expose_value=False,
                 callback=_store_snp_filter,
                 help="Drop SNPs whose minor allele frequency, all populations' alleles pooled, is below this "
                      "(default: 0)"),
    click.option("--maxmaf", type=click.FloatRange(0, 0.5), default=0.5, expose_value=False,
                 callback=_store_snp_filter, help="Drop SNPs whose pooled minor allele frequency is above this "
                                                  "(default: 0.5)"),
    click.option("--keep-snps", type=click.Path(exists=True, dir_okay=False), default=None, expose_value=False,
                 callback=_store_snp_filter,
                 help="Use only the SNP IDs in this file (one per line, or a .snp/.bim file); overrides "
                      "--maxmiss, --minmaf, --maxmaf, --transversions-only and the autosome restriction, as "
                      "admixtools2 keepsnps"),
    click.option("--exclude-snps", type=click.Path(exists=True, dir_okay=False), default=None,
                 expose_value=False, callback=_store_snp_filter,
                 help="Never use the SNP IDs in this file (one per line, or a .snp/.bim file)"),
]


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
        *snp_filter_options,
        *sample_selection_options,
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
        available = set(pm.read_f2_populations(f2_dir))
        missing = [p for p in pops if p not in available]
        if missing:
            raise click.ClickException(f"Not in the f2 blocks in {f2_dir}: {', '.join(missing)}")
        # only these populations' pairs are read
        return pm.read_f2(f2_dir, populations=pops or None)
    return pm.extract_f2(input_prefix, pops, auto_only=not all_chromosomes, block_size=block_size,
                         maxmiss=maxmiss, pseudohaploid=_pseudohaploid(), **snp_filters(),
                         **selected_samples(input_prefix, pops))


def selected_samples(input_prefix: str, pops: list) -> dict:
    """``extract_f2``'s ``inds`` and ``pops`` for the running command's
    sample selection (none without one), POPS checked against the
    populations it leaves."""
    selection = sample_selection()
    if not selection:
        validate_populations(input_prefix, pops)
        return {}
    inds, labels = _core.select_samples(input_prefix, **selection)
    missing = [p for p in pops if p not in set(labels)]
    if missing:
        raise click.ClickException(f"Not among the selected samples' populations: {', '.join(missing)}")
    return {"inds": inds, "pops": labels}


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
    available = pm.read_f2_populations(f2_dir)
    missing = [p for p in pops if p not in set(available)]
    if missing:
        raise click.ClickException(f"Not in the f2 blocks in {f2_dir}: {', '.join(missing)}")
    # only these populations' pairs are read
    return pm.read_f2(f2_dir, populations=pops or None)


def write_table(df, output_format: str, output: str | None, index: bool = False) -> None:
    """Write a result table as TSV or JSON records to OUTPUT or stdout."""
    if output_format == "json":
        text = df.to_json(orient="index" if index else "records", indent=2, double_precision=15) + "\n"
    else:
        text = df.to_csv(sep="\t", index=index, float_format="%.12g")
    if output:
        with open(output, "w") as f:
            f.write(text)
    else:
        click.echo(text, nl=False)


def write_columns(columns: dict, output: str | None, output_format: str = "tsv") -> None:
    """Write equal-length columns as TSV or JSON records, exactly as write_table writes a DataFrame of them.

    For tables built straight from the Rust results, without importing pandas
    (a third of a fast command's time).
    """
    import csv
    import io
    import json

    import numpy as np

    if output_format == "json":
        def value(x):
            if isinstance(x, (float, np.floating)):
                return None if math.isnan(x) else float(x)
            return x.item() if isinstance(x, np.generic) else x

        rows = [dict(zip(columns, map(value, row))) for row in zip(*columns.values())]
        text = json.dumps(rows, indent=2) + "\n"
    else:
        def cell(x) -> str:
            if isinstance(x, (float, np.floating)):
                return "" if math.isnan(x) else f"{x:.12g}"
            return str(x)

        def strings(values) -> list[str]:
            # a column at a time: numpy columns through tolist(), not a cell at a time
            if isinstance(values, np.ndarray) and values.dtype.kind == "f":
                return ["" if v != v else f"{v:.12g}" for v in values.tolist()]
            if isinstance(values, np.ndarray) and values.dtype.kind in "iub":
                return [str(v) for v in values.tolist()]
            return [cell(x) for x in values]

        cols = [strings(v) for v in columns.values()]
        header = list(columns)
        special = ("\t", '"', "\n", "\r")
        if any(c in x for col in [header, *cols] for x in col for c in special):
            # fields csv would quote
            buf = io.StringIO()
            writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(zip(*cols))
            text = buf.getvalue()
        else:
            text = "".join("\t".join(row) + "\n" for row in [header, *zip(*cols)])
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
