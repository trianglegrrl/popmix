"""Commands on samples: PCA, kinship, site-frequency spectra and ADMIXTURE."""

from __future__ import annotations

import click

import popmix as pm

from ._common import (
    _pseudohaploid,
    parse_population_list,
    pseudohaploid_option,
    scan_output_options,
    with_options,
    write_columns,
    write_table,
)
from ._errors import EXPECTED_ERRORS, fail
from ._population_checks import validate_populations


@click.command()
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
@click.option("--shrink", is_flag=True, default=False,
              help="Correct projected samples' shrinkage toward the origin on significant components")
@click.option("--output-prefix", "-o", required=True,
              help="Writes PREFIX.coords.tsv, .eval, .tracy_widom.tsv and .outliers.tsv")
@click.option("--plot", default=None, help="Also plot PC1 against PC2 to this image file (needs popmix[plot])")
def pca(input_prefix: str, pca_pops: list, project_pops: list, components: int, outlier_iterations: int,
        outlier_components: int, outlier_sigma: float, lsqproject: bool, compat: str | None,
        shrink: bool, output_prefix: str, plot: str | None):
    """Principal component analysis as EIGENSOFT smartpca, with projection.

    The PCA is computed on --pca-pops and the samples of --project-pops are
    projected onto it (smartpca poplistname); use --lsqproject for ancient
    samples with missing data.
    """
    try:
        validate_populations(input_prefix, pca_pops + project_pops)
        r = pm.smartpca(input_prefix, pca_pops or None, project_pops, n_components=components,
                        outlier_iterations=outlier_iterations, outlier_components=outlier_components,
                        outlier_sigma=outlier_sigma, lsqproject=lsqproject, compat=compat, shrink=shrink)
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
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
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
@with_options(scan_output_options)
def kinship(input_prefix: str, populations: list, method: str, norm: str, norm_value: float | None,
            window_size: int, compat: str | None, output_format: str, output: str | None):
    """Kinship for every pair of samples: KING-robust or READ-style.

    KING-robust needs diploid genotype calls; for pseudo-haploid ancient
    samples use --method read, which compares them within one group and
    normalizes by the typical (by default median) pair.
    """
    try:
        validate_populations(input_prefix, populations)
        pops = populations or None
        if method == "king":
            from popmix.kinship import king_columns

            write_columns(king_columns(input_prefix, pops), output, output_format)
            return
        else:
            df = pm.read_kinship(input_prefix, pops, norm=norm, norm_value=norm_value, window_size=window_size,
                                 compat=compat)
            click.echo(f"normalization value (P0 of unrelated pairs): {df.attrs['norm_value']:.6g}", err=True)
        write_table(df, output_format, output)
    except click.ClickException:
        raise
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
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
@pseudohaploid_option
@with_options(scan_output_options)
def sfs(input_prefix: str, populations: list, projection: list, folded: bool, joint: bool,
        compat: str | None, all_chromosomes: bool, output_format: str, output: str | None):
    """Site-frequency spectrum, as dadi (projected) or admixtools2 joint_sfs().

    Writes one row per entry: the count in each population and n, the number
    of SNPs (expected, after projection). For one population, also prints S,
    pi, Watterson's theta and Tajima's D to stderr.
    """
    import itertools

    import pandas as pd

    for flag, given in [("--projection", joint and projection), ("--folded", joint and folded),
                        ("--compat", not joint and compat)]:
        if given:
            raise click.UsageError(f"{flag} applies {'without' if joint else 'only with'} --joint")
    try:
        validate_populations(input_prefix, populations)
        if joint:
            df = pm.joint_sfs(input_prefix, populations, compat=compat, auto_only=False)
        else:
            sizes = [int(x) for x in projection] or None
            spec = pm.sfs(input_prefix, populations, sizes, folded=folded, auto_only=not all_chromosomes,
                          pseudohaploid=_pseudohaploid())
            grid = list(itertools.product(*(range(d) for d in spec.data.shape)))
            df = pd.DataFrame(grid, columns=populations)
            df["n"] = spec.data.ravel()
            df["masked"] = spec.mask.ravel()
            if len(populations) == 1 and not folded:
                click.echo(f"S\t{spec.S():.6g}\npi\t{spec.pi():.6g}\n"
                           f"watterson_theta\t{spec.watterson_theta():.6g}\ntajima_d\t{spec.tajima_d():.6g}",
                           err=True)
        write_table(df, output_format, output)
    except click.ClickException:
        raise
    except EXPECTED_ERRORS as e:
        fail(e)


@click.command()
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
    except EXPECTED_ERRORS as e:
        fail(e)
