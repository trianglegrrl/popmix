"""Admixture dating commands: DATES and ROLLOFF."""

from __future__ import annotations

import json
import math
import textwrap

import click

import popmix as pm

from ._common import parse_population_list
from ._errors import EXPECTED_ERRORS, fail
from ._population_checks import validate_populations


@click.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Genotype prefix (with a genetic map)")
@click.option("--admixed", required=True, help="Admixed population")
@click.option("--sources", required=True, callback=parse_population_list,
              help="The two comma-separated source (reference) populations")
@click.option("--binsize", type=float, default=0.0005, help="Distance bin in Morgans (default: 0.0005)")
@click.option("--maxdis", type=float, default=0.05, help="Largest distance in Morgans (default: 0.05)")
@click.option("--runmode", type=click.Choice(["0", "1"]), default="1", help="Weights: 1 pA - pB (default), 0 scaled")
@click.option("--chithresh", type=float, default=6.0,
              help="Leave out pairs in LD in a reference, chi-square at least this (default: 6)")
@click.option("--chromosomes", callback=parse_population_list, default="",
              help="Comma-separated chromosomes to use (default: all autosomes)")
@click.option("--lval", type=float, default=0.5, help="Fit the decay from this distance in cM (default: 0.5)")
@click.option("--hval", type=float, default=100.0, help="Fit the decay up to this distance in cM (default: 100)")
@click.option("--output-prefix", "-o", default=None, help="Also write PREFIX.bins.tsv")
def rolloff(input_prefix: str, admixed: str, sources: list, binsize: float, maxdis: float, runmode: str,
            chithresh: float, chromosomes: list, lval: float, hval: float, output_prefix: str | None):
    """ROLLOFF: date admixture from the decay of weighted LD correlation (ADMIXTOOLS)."""
    try:
        validate_populations(input_prefix, [admixed] + sources)
        r = pm.rolloff(input_prefix, admixed, sources, binsize=binsize, maxdis=maxdis, runmode=int(runmode),
                       chithresh=chithresh, chromosomes=chromosomes or None, lval=lval, hval=hval)
        m, se = r.decay_jackknife
        click.echo(f"decay {r.decay:.3f} generations; jackknife {m:.3f} +/- {se:.3f} "
                   f"({len(r.chromosomes)} chromosomes, {r.n_snps_used} SNPs)")
        if output_prefix:
            r.bins.to_csv(f"{output_prefix}.bins.tsv", sep="\t", index=False, float_format="%.12g")
    except click.ClickException:
        raise
    except EXPECTED_ERRORS as e:
        fail(e)


def _estimate(x: float, se: float, fmt: str) -> str:
    return f"{x:{fmt}}" if math.isnan(se) else f"{x:{fmt}} ± {se:{fmt}}"


def _date_line(r: pm.DatesResult) -> str:
    """The date in years, with the assumptions it rests on (G3)."""
    t = (f"{r.generation_time:g} years/gen, taken as exact" if r.generation_time_se == 0
         else f"{r.generation_time:g} ± {r.generation_time_se:g} years/gen")
    date = _estimate(r.years_bp, r.years_bp_se, ".0f")
    if r.sample_age == 0 and r.sample_age_se == 0:
        text = f"Date: {date} years before the admixed samples lived ({t}; give --sample-age for years BP)"
    else:
        text = f"Date: {date} years BP (sample age {r.sample_age:g} ± {r.sample_age_se:g}, {t})"
    return textwrap.fill(text, width=79, subsequent_indent="      ")


@click.command()
@click.option("--input", "-i", "input_prefix", required=True, help="Input genotype file prefix")
@click.option("--admixed", required=True, help="The admixed population (DATES's test population)")
@click.option("--sources", required=True, callback=parse_population_list,
              help="Its two comma-separated source populations")
@click.option("--bin-width", type=float, default=0.1, help="Distance bin width in cM (default: 0.1)")
@click.option("--max-distance", type=float, default=100.0, help="Maximum genetic distance in cM (default: 100)")
@click.option("--min-distance", type=float, default=0.45,
              help="Fit bins from this distance in cM, DATES lovalfit (default: 0.45)")
@click.option("--generation-time", type=float, default=29.0,
              help="Years per generation, for the date in years (default: 29)")
@click.option("--generation-time-se", type=float, default=0.0,
              help="Standard error of the generation time in years, added to the date's SE (default: 0)")
@click.option("--sample-age", type=float, default=0.0,
              help="Age of the admixed samples in years BP: DATES counts generations before they lived "
                   "(default: 0, present-day)")
@click.option("--sample-age-se", type=float, default=0.0,
              help="Standard error of the sample age in years (default: 0)")
@click.option("--jackknife/--no-jackknife", default=True,
              help="Standard errors by a jackknife over chromosomes (default: on); off is faster, without SEs")
@click.option("--qbin", type=click.IntRange(min=0), default=0,
              help="DATES qbin: 0 visits every SNP pair; N > 0 uses an FFT mesh of N cells "
                   "per bin, faster and approximate (default: 0)")
@click.option("--compat", type=click.Choice(["dates-v753"]), default=None,
              help="Reproduce DATES v753's known quirks instead of correcting them, "
                   "for output identical to DATES (see DEVIATIONS.md)")
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def dates(input_prefix: str, admixed: str, sources: list,
          bin_width: float, max_distance: float, min_distance: float, generation_time: float,
          generation_time_se: float, sample_age: float, sample_age_se: float, jackknife: bool, qbin: int, compat: str | None, output_format: str):
    """Estimate admixture date using DATES algorithm.

    DATES estimates the time since admixture by fitting an exponential decay
    to weighted linkage disequilibrium as a function of genetic distance.
    """
    if len(sources) != 2:
        fail(pm.InputError(f"--sources takes two populations, got {len(sources)}"))
    source1, source2 = sources
    validate_populations(input_prefix, [admixed, source1, source2])

    try:
        click.echo(f"Running DATES: {admixed} = mix of {source1} and {source2}...", err=True)
        result = pm.dates(
            input_prefix,
            admixed,
            sources,
            bin_width_cm=bin_width,
            max_distance_cm=max_distance,
            min_distance_cm=min_distance,
            generation_time=generation_time,
            generation_time_se=generation_time_se,
            sample_age=sample_age,
            sample_age_se=sample_age_se,
            jackknife=jackknife,
            qbin=qbin,
            compat=compat,
        )

        if output_format == "json":
            output = {
                "admixed": admixed,
                "sources": [source1, source2],
                "amplitude": result.fit.amplitude,
                "generations": result.generations,
                "generations_se": result.generations_se,
                "years_bp": result.years_bp,
                "years_bp_se": result.years_bp_se,
                "generation_time": result.generation_time,
                "generation_time_se": result.generation_time_se,
                "sample_age": result.sample_age,
                "sample_age_se": result.sample_age_se,
                "n_snps": result.n_snps,
                "source_alleles": list(result.source_alleles),
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
            click.echo(f"\nDATES Results: {admixed} = mix of {source1} and {source2}")
            click.echo("=" * 60)
            click.echo(f"Amplitude (fit): {result.fit.amplitude:.6f}")
            click.echo(f"Estimated generations: {_estimate(result.generations, result.generations_se, '.1f')}")
            click.echo(_date_line(result))
            click.echo(f"Fit R²: {result.fit.r_squared:.4f} (how well one exponential fits; it does not show a "
                       "single pulse)")
            click.echo(f"SNPs used: {result.n_snps}")
            click.echo(f"LD bins: {len(result.ld_bins)}")
            click.echo()
            click.echo("LD Decay (first 10 bins):")
            click.echo("  Dist(cM)  Pairs      Correlation")
            for b in result.ld_bins[:10]:
                click.echo(f"  {b.distance_cm:6.1f}    {b.n_pairs:9d}  {b.correlation:+.6f}")

    except EXPECTED_ERRORS as e:
        fail(e)
