//! Some of a dataset's samples and SNPs, read a block of samples at a time
//! for the sample-major writers (PED, TGENO), from a file read in place
//! (TASKS.md P4): `convertf --populations ... -of tgeno` of AADR holds a
//! block of the kept samples, not the dataset.

use rayon::prelude::*;

use crate::streaming::SnpMajor;
use crate::tgeno::TgenoReader;
use crate::writers::SampleGenotypes;

/// A file read in place.
pub enum Source<'a> {
    SnpMajor(&'a SnpMajor),
    Tgeno(&'a TgenoReader),
}

/// The samples `samples` and SNPs `snps` (indices in the file, in output
/// order) of `source`.
pub struct SubsetSamples<'a> {
    pub source: Source<'a>,
    pub samples: Vec<usize>,
    pub snps: Vec<usize>,
}

impl SampleGenotypes for SubsetSamples<'_> {
    fn fill_samples(&self, samples: std::ops::Range<usize>, n_snps: usize, out: &mut [u8]) {
        let of = &self.samples[samples];
        let snps = &self.snps;
        match self.source {
            Source::SnpMajor(file) => file.fill_samples_of(of, n_snps, out, &|j| (snps[j], false)),
            Source::Tgeno(reader) => {
                out.par_chunks_mut(n_snps.max(1))
                    .zip(of)
                    .for_each(|(row, &s)| {
                        let record = reader.record(s);
                        for (o, &snp) in row.iter_mut().zip(snps) {
                            *o = (record[snp / 4] >> (6 - 2 * (snp % 4))) & 0b11;
                        }
                    });
                // the kept samples are in file order: those before are done
                if let Some(&last) = of.last() {
                    if of.windows(2).all(|w| w[0] < w[1]) {
                        reader.release_before(last + 1);
                    }
                }
            }
        }
    }
}
