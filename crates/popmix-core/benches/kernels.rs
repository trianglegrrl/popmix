//! Benchmarks of popmix's hot kernels on synthetic data (TASKS.md F12):
//! `cargo bench -p popmix-core --bench kernels`.
//!
//! The data are random genotypes with some missing calls; sizes are small
//! enough to run in seconds and large enough for the kernels to dominate.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array2;

use popmix_core::{
    admixture, king_robust, smartpca, AdmixtureOptions, AlleleFrequencies, BitGenotypes,
    BlockDefinition, Chromosome, F2Blocks, F2Options, GenotypeMatrix, MatrixRows, PackedGenotypes,
    PcaOptions, PopulationId, Sample, SampleRole, Snp,
};

/// Genotypes of `n` samples in four populations at `m` SNPs, 5% missing.
fn genotypes(n: usize, m: usize) -> GenotypeMatrix {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let samples = (0..n)
        .map(|i| Sample::new(format!("s{i}"), format!("P{}", i % 4)))
        .collect();
    let snps = (0..m)
        .map(|j| {
            Snp::new(
                format!("rs{j}"),
                Chromosome::Autosome(1 + (j * 22 / m) as u8),
                j as u32,
                j as f64 * 1e-5,
                'A',
                'G',
            )
        })
        .collect();
    let data = Array2::from_shape_fn((m, n), |_| match next() % 20 {
        0 => 3,
        x => (x % 3) as u8,
    });
    GenotypeMatrix::new(data, samples, snps).unwrap()
}

fn rows(g: &GenotypeMatrix) -> MatrixRows<'_> {
    MatrixRows::new(g, (0..g.n_snps()).collect(), (0..g.n_samples()).collect())
}

fn pca(c: &mut Criterion) {
    let mut group = c.benchmark_group("pca");
    group.sample_size(10);
    for n in [200, 500] {
        let g = genotypes(n, 20_000);
        let packed = PackedGenotypes::read(&mut rows(&g), n, g.n_snps()).unwrap();
        let roles = vec![SampleRole::Pca; n];
        let opts = {
            let mut o = PcaOptions::default();
            o.outlier_iterations = 0;
            o
        };
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| smartpca(&packed, &roles, &opts).unwrap())
        });
    }
    group.finish();
}

fn king(c: &mut Criterion) {
    let mut group = c.benchmark_group("king");
    group.sample_size(10);
    for n in [100, 300] {
        let g = genotypes(n, 50_000);
        let bits = BitGenotypes::read(&mut rows(&g), n, g.n_snps()).unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| king_robust(&bits))
        });
    }
    group.finish();
}

fn f2_blocks(c: &mut Criterion) {
    let mut group = c.benchmark_group("f2_blocks");
    group.sample_size(10);
    let g = genotypes(400, 50_000);
    let af = AlleleFrequencies::from_genotypes(&g);
    let opts = {
        let mut o = F2Options::default();
        o.max_missing = 1.0;
        o.poly_only = true;
        o
    };
    let def = BlockDefinition::new(0.005);
    group.bench_function("4 populations, maxmiss 1", |b| {
        b.iter(|| F2Blocks::from_allele_frequencies_with_options(&af, &def, &opts))
    });
    // many populations: one per pair of samples
    let k = 100;
    let wide = AlleleFrequencies::new(
        Array2::from_shape_fn((af.n_snps(), k), |(s, p)| {
            (af.frequencies()[[s, p % 4]] + p as f64 * 1e-3).min(1.0)
        }),
        Array2::from_shape_fn((af.n_snps(), k), |(s, p)| af.counts()[[s, p % 4]]),
        (0..k).map(|p| PopulationId::new(format!("W{p}"))).collect(),
        af.snps().to_vec(),
    )
    .unwrap();

    group.bench_function("100 populations, maxmiss 1", |b| {
        b.iter(|| F2Blocks::from_allele_frequencies_with_options(&wide, &def, &opts))
    });
    group.finish();
}

fn admixture_fit(c: &mut Criterion) {
    let mut group = c.benchmark_group("admixture");
    group.sample_size(10);
    let g = genotypes(60, 5_000);
    let packed = PackedGenotypes::read(&mut rows(&g), 60, g.n_snps()).unwrap();
    let opts = {
        let mut o = AdmixtureOptions::default();
        o.max_iterations = 5;
        o
    };
    group.bench_function("K=3, 5 iterations", |b| {
        b.iter(|| admixture(&packed, 3, &opts).unwrap())
    });
    group.finish();
}

criterion_group!(kernels, pca, king, f2_blocks, admixture_fit);
criterion_main!(kernels);
