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
"""

import json
import sys
from typing import Optional, Tuple

import click

import popmix as pm


def parse_population_list(ctx, param, value: str) -> list:
    """Parse comma-separated population list."""
    if not value:
        return []
    return [p.strip() for p in value.split(",") if p.strip()]


def parse_snp_range(ctx, param, value: Optional[str]) -> Optional[Tuple[int, int]]:
    """Parse SNP range in format 'start:end'."""
    if not value:
        return None
    try:
        parts = value.split(":")
        if len(parts) != 2:
            raise click.BadParameter("SNP range must be in format 'start:end'")
        return (int(parts[0]), int(parts[1]))
    except ValueError:
        raise click.BadParameter("SNP range must be integers in format 'start:end'")


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
            click.echo(f"\nTip: Use --populations to show all populations with counts")
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("prefix")
@click.option("--filter", "-f", "pattern", help="Filter populations (supports wildcards like French*)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def populations(prefix: str, pattern: Optional[str], output_format: str):
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
                    click.echo(f'Tip: Use --filter "pattern*" to search')
                    click.echo(f'Tip: Specify populations with: --populations "{",".join(example_pops)}"')

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command("extract-f2")
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--output", "-o", "output_dir", required=True, help="Output directory for F2 blocks")
@click.option("--snp-range", callback=parse_snp_range, help="SNP range to process (format: start:end)")
@click.option("--block-size", type=float, default=0.05, help="Block size in Morgans (default: 0.05)")
@click.option("--populations", callback=parse_population_list, help="Comma-separated list of populations to include")
@click.option("--auto-only/--all-chromosomes", default=True, help="Use only autosomal SNPs (default: auto-only)")
def extract_f2(input_prefix: str, output_dir: str, snp_range: Optional[Tuple[int, int]],
               block_size: float, populations: list, auto_only: bool):
    """Extract F2 blocks from genotype data.

    Computes pairwise F2 statistics for all population pairs and saves them
    for use in downstream analyses (f3, f4, qpAdm, etc.).
    """
    import os
    try:
        # Use efficient streaming extraction when populations are specified
        # and no SNP range filter is needed
        if populations and not snp_range:
            click.echo(f"Streaming F2 extraction from {input_prefix}...")
            click.echo(f"Populations: {', '.join(populations)}")
            f2_blocks = pm.extract_f2(
                input_prefix,
                populations,
                auto_only=auto_only,
                block_size=block_size
            )
            click.echo(f"Computed {f2_blocks.n_blocks} blocks for {f2_blocks.n_populations} populations")
        else:
            # Fall back to loading genotypes when SNP range is specified
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
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--pop1", required=True, help="First population")
@click.option("--pop2", required=True, help="Second population")
@click.option("--snp-range", callback=parse_snp_range, help="SNP range to process (format: start:end)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def f2(input_prefix: str, pop1: str, pop2: str, snp_range: Optional[Tuple[int, int]], output_format: str):
    """Compute f2 statistic between two populations.

    f2(A, B) measures genetic drift between populations A and B.
    """
    # Validate populations early
    validate_populations(input_prefix, [pop1, pop2])

    try:
        geno = pm.read_genotypes(input_prefix, snp_range=snp_range)
        # Compute frequencies only for the populations we need (much faster)
        af = geno.allele_frequencies(populations=[pop1, pop2])
        f2_blocks = pm.f2_from_frequencies(af)

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


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--target", required=True, help="Target population (potentially admixed)")
@click.option("--source1", required=True, help="First source population")
@click.option("--source2", required=True, help="Second source population")
@click.option("--snp-range", callback=parse_snp_range, help="SNP range to process (format: start:end)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def f3(input_prefix: str, target: str, source1: str, source2: str,
       snp_range: Optional[Tuple[int, int]], output_format: str):
    """Compute f3 statistic (admixture test).

    f3(A; B, C) tests whether A is admixed between B and C.
    A negative value indicates admixture.
    """
    # Validate populations early
    validate_populations(input_prefix, [target, source1, source2])

    try:
        geno = pm.read_genotypes(input_prefix, snp_range=snp_range)
        # Compute frequencies only for the populations we need (much faster)
        af = geno.allele_frequencies(populations=[target, source1, source2])
        f2_blocks = pm.f2_from_frequencies(af)

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
                click.echo(f"  => Significant admixture signal")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--pop1", required=True, help="Population 1")
@click.option("--pop2", required=True, help="Population 2")
@click.option("--pop3", required=True, help="Population 3")
@click.option("--pop4", required=True, help="Population 4")
@click.option("--snp-range", callback=parse_snp_range, help="SNP range to process (format: start:end)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def f4(input_prefix: str, pop1: str, pop2: str, pop3: str, pop4: str,
       snp_range: Optional[Tuple[int, int]], output_format: str):
    """Compute f4 statistic (D-statistic / ABBA-BABA test).

    f4(A, B; C, D) tests for gene flow between (A,B) and (C,D) branches.
    """
    # Validate populations early
    validate_populations(input_prefix, [pop1, pop2, pop3, pop4])

    try:
        geno = pm.read_genotypes(input_prefix, snp_range=snp_range)
        # Compute frequencies only for the populations we need (much faster)
        af = geno.allele_frequencies(populations=[pop1, pop2, pop3, pop4])
        f2_blocks = pm.f2_from_frequencies(af)

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
                click.echo(f"  => Significant deviation from tree-like evolution")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--target", required=True, help="Target population to model")
@click.option("--left", required=True, callback=parse_population_list, help="Source populations (comma-separated)")
@click.option("--right", required=True, callback=parse_population_list, help="Outgroup populations (comma-separated)")
@click.option("--snp-range", callback=parse_snp_range, help="SNP range to process (format: start:end)")
@click.option("--constrained/--unconstrained", default=False, help="Constrain weights to [0, 1]")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpadm(input_prefix: str, target: str, left: list, right: list,
          snp_range: Optional[Tuple[int, int]], constrained: bool, output_format: str):
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

    # Validate populations early
    all_pops = [target] + list(left) + list(right)
    validate_populations(input_prefix, all_pops)

    try:
        geno = pm.read_genotypes(input_prefix, snp_range=snp_range)
        # Compute frequencies only for the populations we need (much faster)
        af = geno.allele_frequencies(populations=all_pops)
        f2_blocks = pm.f2_from_frequencies(af)

        result = pm.qpadm(f2_blocks, target, left, right, constrained=constrained)

        if output_format == "json":
            output = {
                "target": result.target,
                "sources": result.sources,
                "weights": dict(zip(result.sources, result.weights)),
                "chi_sq": result.chi_sq,
                "df": result.df,
                "p_value": result.p_value,
                "feasible": result.is_feasible(),
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("target\tsource\tweight\tchi_sq\tdf\tp_value")
            for source, weight in zip(result.sources, result.weights):
                click.echo(f"{target}\t{source}\t{weight:.4f}\t{result.chi_sq:.4f}\t{result.df}\t{result.p_value:.4f}")
        else:
            click.echo(f"qpAdm: {target} = mix of {', '.join(left)}")
            click.echo(f"Outgroups: {', '.join(right)}")
            click.echo()
            click.echo("Admixture weights:")
            for source, weight in zip(result.sources, result.weights):
                click.echo(f"  {source}: {weight:.4f} ({weight*100:.1f}%)")
            click.echo()
            click.echo(f"Model fit: chi-sq={result.chi_sq:.4f}, df={result.df}, p={result.p_value:.4f}")
            click.echo(f"Feasible: {result.is_feasible()}")
            if result.is_good_fit():
                click.echo("=> Good model fit (p > 0.05)")
            else:
                click.echo("=> Poor model fit (p <= 0.05)")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--left", required=True, callback=parse_population_list, help="Left populations (comma-separated)")
@click.option("--right", required=True, callback=parse_population_list, help="Right populations (comma-separated)")
@click.option("--snp-range", callback=parse_snp_range, help="SNP range to process (format: start:end)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpwave(input_prefix: str, left: list, right: list,
           snp_range: Optional[Tuple[int, int]], output_format: str):
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

    # Validate populations early
    all_pops = list(left) + list(right)
    validate_populations(input_prefix, all_pops)

    try:
        geno = pm.read_genotypes(input_prefix, snp_range=snp_range)
        # Compute frequencies only for the populations we need (much faster)
        af = geno.allele_frequencies(populations=all_pops)
        f2_blocks = pm.f2_from_frequencies(af)

        result = pm.qpwave(f2_blocks, left, right)

        if output_format == "json":
            output = {
                "left": result.left_pops,
                "right": result.right_pops,
                "estimated_rank": result.estimated_rank,
                "rank_tests": [
                    {"rank": t.rank, "chi_sq": t.chi_sq, "df": t.df, "p_value": t.p_value}
                    for t in result.rank_tests
                ],
                "is_clade": result.is_clade(),
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("rank\tchi_sq\tdf\tp_value")
            for t in result.rank_tests:
                click.echo(f"{t.rank}\t{t.chi_sq:.4f}\t{t.df}\t{t.p_value:.6f}")
        else:
            click.echo(f"qpWave: {', '.join(left)}")
            click.echo(f"Outgroups: {', '.join(right)}")
            click.echo()
            click.echo("Rank tests:")
            click.echo("  Rank  Chi-sq     df    P-value")
            for t in result.rank_tests:
                sig = "*" if t.p_value < 0.05 else " "
                click.echo(f"  {t.rank:4d}  {t.chi_sq:8.4f}  {t.df:4d}  {t.p_value:.6f} {sig}")
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
def convertf(input_prefix: str, output_prefix: str, input_format: Optional[str], output_format: str):
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
          input_format_a: Optional[str], input_format_b: Optional[str]):
    """Merge two genotype datasets.

    Combines samples from two genotype files, keeping only SNPs that
    exist in both datasets. Alleles are matched and genotypes are
    recoded when REF/ALT are swapped between files.

    Supports VCF input (pass full path with .vcf or .vcf.gz extension).
    VCF output is not supported.

    Based on PLINK 2.0 merge algorithm by Shaun Purcell & Christopher Chang.
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

        if result.get('flipped_snps', 0) > 0:
            click.echo(f"  Flipped SNPs (recoded): {result['flipped_snps']}")
        if result.get('strand_ambiguous_snps', 0) > 0:
            click.echo(f"  Strand-ambiguous SNPs (skipped): {result['strand_ambiguous_snps']}")

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
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def dates(input_prefix: str, test_pop: str, source1: str, source2: str,
          bin_width: float, max_distance: float, output_format: str):
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


def main():
    """Entry point for the popmix CLI."""
    cli()


if __name__ == "__main__":
    main()
