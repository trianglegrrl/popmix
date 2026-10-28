"""Commands about datasets: what they hold, converting and merging them, and the agent skill."""

from __future__ import annotations

import click

import popmix as pm

from ._common import (
    _pairs,
    _sample_ids,
    _sample_selection_options,
    _snp_ids,
    parse_population_list,
    with_options,
)
from ._errors import EXPECTED_ERRORS, fail


def _parse_region(ctx, param, value: str | None) -> tuple[str, int, int] | None:
    """CHROM:START-END as (chromosome, start, end)."""
    if not value:
        return None
    try:
        chrom, span = value.rsplit(":", 1)
        start, end = (int(x.replace(",", "")) for x in span.split("-"))
    except ValueError:
        raise click.BadParameter("give the region as CHROM:START-END, e.g. 2:135000000-137000000") from None
    if start > end:
        raise click.BadParameter(f"the region starts after it ends ({start} > {end})")
    return chrom, start, end


@click.command()
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
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
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

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--par", "-p", "par_file", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Run a convertf parameter file (genotypename, snpname, indivname, outputformat, "
                   "genotypeoutname, ..., poplistname, badsnpname, chrom, lopos, hipos) instead of the options")
@click.option("--input", "-i", "input_prefix", default=None, help="Input genotype file prefix (or full path for VCF)")
@click.option("--output", "-o", "output_prefix", default=None, help="Output genotype file prefix")
@click.option("--input-format", "-if", "input_format", default=None,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap", "tgeno",
                                "vcf"]),
              help="Input format (auto-detected if not specified)")
@click.option("--output-format", "-of", "output_format", default=None,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap",
                                "tgeno", "vcf"]),
              help="Output format")
@with_options(_sample_selection_options(exposed=True))
@click.option("--populations", callback=parse_population_list, default="",
              help="Comma-separated populations to keep (after --relabel and --relabel-pops)")
@click.option("--keep-snps", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Keep only the SNP IDs in this file (one per line, or a .snp/.bim file)")
@click.option("--exclude-snps", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Leave out the SNP IDs in this file (one per line, or a .snp/.bim file)")
@click.option("--chromosomes", callback=parse_population_list, default="",
              help="Comma-separated chromosomes to keep (1-22, X, Y, MT)")
@click.option("--region", callback=_parse_region, default=None,
              help="Keep only SNPs in CHROM:START-END (base pairs, both ends included)")
def convertf(par_file: str | None, input_prefix: str | None, output_prefix: str | None, input_format: str | None,
             output_format: str | None,
             keep_inds: str | None, exclude_inds: str | None, relabel: str | None, relabel_pops: str | None,
             populations: list, keep_snps: str | None, exclude_snps: str | None, chromosomes: list,
             region: tuple[str, int, int] | None):
    """Convert between genotype file formats.

    Reads and writes eight formats:

    \b
    - vcf: VCF (.vcf, .vcf.gz), given by its full path; written as PREFIX.vcf
    - ancestrymap: sparse text (.ancestrymapgeno, .snp, .ind)
    - eigenstrat: dense text (.geno, .snp, .ind)
    - packedancestrymap: binary EIGENSTRAT (.geno, .snp, .ind)
    - tgeno: packed genotypes stored sample by sample (.geno, .snp, .ind; AADR v66 on)
    - packedped/plink: binary PLINK (.bed, .bim, .fam)
    - ped: text PLINK (.ped, .map, .pedind)

    Samples can be kept or left out by ID, given new populations and
    chosen by population; SNPs by ID, chromosome and region.
    """
    if par_file is not None:
        given = [k for k, v in [("--input", input_prefix), ("--output", output_prefix),
                                ("--output-format", output_format)] if v is not None]
        if given:
            raise click.UsageError(f"{', '.join(given)}: give either --par or the options, not both")
        try:
            from popmix.convertf_par import convertf_par

            result = convertf_par(par_file)
            click.echo(f"Converted {result['n_snps']} SNPs x {result['n_samples']} samples")
        except EXPECTED_ERRORS as e:
            fail(e)
        return
    for flag, value in [("--input", input_prefix), ("--output", output_prefix), ("--output-format", output_format)]:
        if value is None:
            raise click.UsageError(f"Missing option {flag} (or give --par)")
    try:
        # Detect input format if not specified
        if input_format is None:
            detected = pm.detect_format(input_prefix)
            click.echo(f"Auto-detected input format: {detected}")
        else:
            click.echo(f"Input format: {input_format}")

        click.echo(f"Output format: {output_format}")
        click.echo(f"Converting {input_prefix} -> {output_prefix}...")

        subset = {
            "keep_inds": keep_inds and _sample_ids(keep_inds),
            "exclude_inds": exclude_inds and _sample_ids(exclude_inds),
            "relabel": relabel and _pairs(relabel),
            "relabel_pops": relabel_pops and _pairs(relabel_pops),
            "populations": populations or None,
            "keep_snps": keep_snps and _snp_ids(keep_snps),
            "exclude_snps": exclude_snps and _snp_ids(exclude_snps),
            "chromosomes": chromosomes or None,
            "region": region,
        }
        result = pm.convertf(
            input_prefix=input_prefix,
            output_prefix=output_prefix,
            output_format=output_format,
            input_format=input_format,
            **{k: v for k, v in subset.items() if v is not None},
        )

        click.echo(f"Converted {result['n_snps']} SNPs x {result['n_samples']} samples")
        click.echo(f"Output files written to {output_prefix}")

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
@click.option("--file-a", "-a", required=True, help="First input genotype file (prefix or full VCF path)")
@click.option("--file-b", "-b", required=True, help="Second input genotype file (prefix or full VCF path)")
@click.option("--output", "-o", required=True, help="Output genotype file prefix")
@click.option("--output-format", "-of", required=True,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap",
                                "tgeno", "vcf"]),
              help="Output format")
@click.option("--input-format-a", "-ifa", default=None,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap", "tgeno",
                                "vcf"]),
              help="Format for file A (auto-detected if not specified)")
@click.option("--input-format-b", "-ifb", default=None,
              type=click.Choice(["ancestrymap", "eigenstrat", "ped", "packedped", "plink", "packedancestrymap", "tgeno",
                                "vcf"]),
              help="Format for file B (auto-detected if not specified)")
@click.option("--strandcheck/--no-strandcheck", default=True,
              help="Reconcile strands: complemented alleles are the other strand, A/T and C/G SNPs are left out "
                   "(default). --no-strandcheck when both files are on one strand (two releases of one panel): "
                   "A/T and C/G SNPs are kept, complemented alleles are a mismatch (mergeit strandcheck: NO)")
def merge(file_a: str, file_b: str, output: str, output_format: str,
          input_format_a: str | None, input_format_b: str | None, strandcheck: bool):
    """Merge two genotype datasets.

    Combines samples from two genotype files, keeping only SNPs that
    exist in both datasets, paired by position (and by name where several
    share one). Alleles are reconciled as AdmixTools mergeit does: genotypes
    are recoded when the alleles are swapped, a SNP on the other strand is
    complemented, the unknown allele of a monomorphic SNP is filled in, and
    A/T, C/G and third-allele SNPs are left out and counted.

    Supports VCF input (pass full path with .vcf or .vcf.gz extension).
    VCF output is written as OUTPUT.vcf (OUTPUT itself if it ends in .vcf
    or .vcf.gz, compressed for .gz).
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
            strandcheck=strandcheck,
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

    except EXPECTED_ERRORS as e:
        fail(e)


SKILL_DEST = "~/.claude/skills/popmix"


def _is_popmix_skill(path) -> bool:
    """Whether PATH holds a copy of this skill (so --force may replace it)."""
    skill_md = path / "SKILL.md"
    try:
        return skill_md.is_file() and "\nname: popmix\n" in skill_md.read_text(errors="replace")
    except OSError:
        return False


def _copy_resource_tree(source, dest) -> int:
    """Copy a package resource directory to DEST; the number of files copied."""
    dest.mkdir(parents=True)
    n = 0
    for item in sorted(source.iterdir(), key=lambda p: p.name):
        if item.name == "__pycache__":
            continue
        if item.is_dir():
            n += _copy_resource_tree(item, dest / item.name)
        else:
            (dest / item.name).write_bytes(item.read_bytes())
            n += 1
    return n


@click.command("install-skill")
@click.option("--dest", default=SKILL_DEST,
              help=f"Directory to install the skill into (default: {SKILL_DEST}, Claude Code's "
                   "personal skills)")
@click.option("--force", is_flag=True, default=False,
              help="Replace an earlier copy of the popmix skill at --dest")
def install_skill(dest: str, force: bool):
    """Install the popmix agent skill (SKILL.md and one file per command) for Claude Code.

    The skill ships inside the popmix package and documents this version's
    commands, so reinstall it (with --force) after upgrading popmix.
    """
    import shutil
    from importlib import resources
    from pathlib import Path

    try:
        target = Path(dest).expanduser()
        if target.exists() or target.is_symlink():
            if not target.is_dir():
                raise ValueError(f"{target} exists and is not a directory")
            if any(target.iterdir()):
                if not _is_popmix_skill(target):
                    raise ValueError(f"{target} is not empty and not a popmix skill; choose another --dest")
                if not force:
                    raise ValueError(f"{target} already holds a popmix skill; use --force to replace it")
                shutil.rmtree(target)
            else:
                target.rmdir()
        n = _copy_resource_tree(resources.files("popmix") / "skill", target)
        click.echo(f"Installed the popmix {pm.__version__} skill ({n} files) to {target}")
    except EXPECTED_ERRORS as e:
        fail(e)
