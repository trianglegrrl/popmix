//! Genetic map block definitions for jackknife resampling.
//!
//! SNPs are grouped into blocks based on genetic map position (in Morgans).
//! The default block size is 0.05 Morgans (5 centiMorgans). As admixtools2
//! `get_block_lengths`: a block size of 100 or more is in base pairs, and
//! without a genetic map (fewer than two distinct map positions, as in a
//! VCF) blocks are 2 Mb of base-pair position.

use crate::genotypes::Snp;

/// Default block size in Morgans (5 cM).
pub const DEFAULT_BLOCK_SIZE_MORGAN: f64 = 0.05;

/// Block size in base pairs when there is no genetic map (admixtools2).
pub const NO_MAP_BLOCK_SIZE_BP: f64 = 2e6;

/// Definition of how to partition SNPs into blocks.
#[derive(Debug, Clone)]
pub struct BlockDefinition {
    /// Block size in Morgans.
    pub block_size: f64,
}

impl Default for BlockDefinition {
    fn default() -> Self {
        BlockDefinition {
            block_size: DEFAULT_BLOCK_SIZE_MORGAN,
        }
    }
}

impl BlockDefinition {
    /// Create a new block definition with custom block size.
    pub fn new(block_size_morgan: f64) -> Self {
        BlockDefinition {
            block_size: block_size_morgan,
        }
    }

    /// The coordinate blocks of SNPs with this `evidence` are measured in.
    pub fn coordinate(&self, evidence: &CoordinateEvidence) -> BlockCoordinate {
        if self.block_size >= 100.0 {
            BlockCoordinate::BasePairs(self.block_size)
        } else if evidence.morgans_differ == 0 && evidence.bp_differ >= 2 {
            BlockCoordinate::BasePairs(NO_MAP_BLOCK_SIZE_BP)
        } else {
            BlockCoordinate::Morgans(self.block_size)
        }
    }

    /// Assign SNPs to blocks based on their genetic map positions.
    ///
    /// Returns a BlockAssignment that maps each SNP index to a block index.
    pub fn assign_blocks(&self, snps: &[Snp]) -> BlockAssignment {
        let mut evidence = CoordinateEvidence::default();
        evidence.observe(snps);
        self.assign_blocks_in(snps, self.coordinate(&evidence))
    }

    /// [`assign_blocks`](Self::assign_blocks) in the coordinate given. A
    /// block never spans two chromosomes, so the blocks of a run of whole
    /// chromosomes are those of the run on its own.
    pub fn assign_blocks_in(&self, snps: &[Snp], coordinate: BlockCoordinate) -> BlockAssignment {
        if snps.is_empty() {
            return BlockAssignment {
                snp_to_block: Vec::new(),
                block_starts: Vec::new(),
                block_lengths: Vec::new(),
                n_blocks: 0,
            };
        }
        let bp = |s: &Snp| s.position_bp as f64;
        let morgans = |s: &Snp| s.position_morgans;
        let (size, coord): (f64, &dyn Fn(&Snp) -> f64) = match coordinate {
            BlockCoordinate::BasePairs(size) => (size, &bp),
            BlockCoordinate::Morgans(size) => (size, &morgans),
        };

        let mut snp_to_block = Vec::with_capacity(snps.len());
        let mut block_starts = Vec::new();
        let mut block_lengths = Vec::new();

        let mut current_block = 0;
        let mut current_block_start = 0;
        let mut current_chr = snps[0].chromosome.clone();
        let mut current_block_start_pos = coord(&snps[0]);
        let mut previous_pos = current_block_start_pos;

        for (snp_idx, snp) in snps.iter().enumerate() {
            // Start a new block if:
            // 1. Chromosome changed
            // 2. Exceeded block size in genetic distance (block_size is in Morgans)
            // 3. The position went back (unsorted input): a block is a run of
            //    neighbouring SNPs
            let pos = coord(snp);
            let new_block = snp.chromosome != current_chr
                || (pos - current_block_start_pos) >= size
                || pos < previous_pos;
            previous_pos = pos;

            if new_block && snp_idx > 0 {
                // Record the previous block
                block_starts.push(current_block_start);
                block_lengths.push(snp_idx - current_block_start);

                // Start new block
                current_block += 1;
                current_block_start = snp_idx;
                current_chr = snp.chromosome.clone();
                current_block_start_pos = coord(snp);
            }

            snp_to_block.push(current_block);
        }

        // Record the final block
        block_starts.push(current_block_start);
        block_lengths.push(snps.len() - current_block_start);

        BlockAssignment {
            snp_to_block,
            block_starts,
            block_lengths,
            n_blocks: current_block + 1,
        }
    }
}

/// The coordinate blocks are measured in, with their size in it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BlockCoordinate {
    BasePairs(f64),
    Morgans(f64),
}

/// What decides the coordinate (admixtools2: base pairs when there is no
/// genetic map, fewer than two distinct map positions), gathered over runs of
/// SNPs in order: how many SNPs differ from the first one's map and base-pair
/// positions, up to two.
#[derive(Debug, Clone, Default)]
pub struct CoordinateEvidence {
    first: Option<(f64, u32)>,
    morgans_differ: usize,
    bp_differ: usize,
}

impl CoordinateEvidence {
    /// Take the next run of SNPs into account.
    pub fn observe(&mut self, snps: &[Snp]) {
        let Some(snp) = snps.first() else { return };
        let (morgans, bp) = *self
            .first
            .get_or_insert((snp.position_morgans, snp.position_bp));
        let more = |have: usize, differs: &dyn Fn(&Snp) -> bool| {
            have + snps.iter().filter(|s| differs(s)).take(2 - have).count()
        };
        if self.morgans_differ < 2 {
            self.morgans_differ = more(self.morgans_differ, &|s| s.position_morgans != morgans);
        }
        if self.bp_differ < 2 {
            self.bp_differ = more(self.bp_differ, &|s| s.position_bp != bp);
        }
    }
}

/// Assignment of SNPs to blocks.
#[derive(Debug, Clone)]
pub struct BlockAssignment {
    /// For each SNP index, which block it belongs to.
    pub snp_to_block: Vec<usize>,
    /// Starting SNP index for each block.
    pub block_starts: Vec<usize>,
    /// Number of SNPs in each block.
    pub block_lengths: Vec<usize>,
    /// Total number of blocks.
    pub n_blocks: usize,
}

impl BlockAssignment {
    /// Get the block index for a given SNP.
    #[inline]
    pub fn block_of(&self, snp_idx: usize) -> usize {
        self.snp_to_block[snp_idx]
    }

    /// Get the SNP range for a given block.
    pub fn block_range(&self, block_idx: usize) -> std::ops::Range<usize> {
        let start = self.block_starts[block_idx];
        let end = start + self.block_lengths[block_idx];
        start..end
    }

    /// Get the number of SNPs in a block.
    #[inline]
    pub fn block_len(&self, block_idx: usize) -> usize {
        self.block_lengths[block_idx]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::Chromosome;

    fn make_snps_same_chr(n: usize, cm_spacing: f64) -> Vec<Snp> {
        // Convert cM to Morgans for internal storage
        let morgans_spacing = cm_spacing / 100.0;
        (0..n)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    i as f64 * morgans_spacing,
                    'A',
                    'G',
                )
            })
            .collect()
    }

    #[test]
    fn test_large_block_sizes_are_base_pairs() {
        // SNPs every 1000 bp: 2500 bp blocks hold 3 SNPs
        let snps = make_snps_same_chr(9, 1.0);
        let blocks = BlockDefinition::new(2500.0).assign_blocks(&snps);
        assert_eq!(blocks.block_lengths, vec![3, 3, 3]);
    }

    #[test]
    fn test_no_genetic_map_uses_2mb_blocks() {
        let snps: Vec<Snp> = (0..10)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    i * 500_000,
                    0.0,
                    'A',
                    'G',
                )
            })
            .collect();
        let blocks = BlockDefinition::default().assign_blocks(&snps);
        assert_eq!(blocks.block_lengths, vec![4, 4, 2]);
    }

    #[test]
    fn a_position_going_back_starts_a_block() {
        // cM 0, 1, 2, then back to 0.5: the last SNP is not next to the others
        let mut snps = make_snps_same_chr(4, 1.0);
        snps[3].position_morgans = 0.005;
        let blocks = BlockDefinition::default().assign_blocks(&snps);
        assert_eq!(blocks.block_lengths, vec![3, 1]);
    }

    #[test]
    fn test_single_block_small_region() {
        // 10 SNPs at 1 cM spacing in 5 cM blocks: 0-4 cM in block 0, 5-9 cM
        // in block 1
        let snps = make_snps_same_chr(10, 1.0);
        let blocks = BlockDefinition::default().assign_blocks(&snps);

        assert_eq!(blocks.n_blocks, 2);
        assert_eq!(blocks.block_of(0), 0);
        assert_eq!(blocks.block_of(4), 0);
        assert_eq!(blocks.block_of(5), 1);
        assert_eq!(blocks.block_of(9), 1);
    }

    #[test]
    fn test_multiple_blocks() {
        // 100 SNPs at 1 cM spacing = 100 cM total
        // With 5 cM blocks, should have ~20 blocks (may vary slightly due to FP precision)
        let snps = make_snps_same_chr(100, 1.0);
        let blocks = BlockDefinition::default().assign_blocks(&snps);

        // Allow some tolerance for floating point precision at boundaries
        assert!(
            blocks.n_blocks >= 18 && blocks.n_blocks <= 22,
            "Expected ~20 blocks, got {}",
            blocks.n_blocks
        );

        // Check that all SNPs are assigned
        assert_eq!(blocks.snp_to_block.len(), 100);
    }

    #[test]
    fn test_chromosome_boundary_creates_new_block() {
        // Create SNPs on two chromosomes
        let mut snps = make_snps_same_chr(5, 1.0);
        for i in 0..5 {
            snps.push(Snp::new(
                format!("rs{}", 100 + i),
                Chromosome::Autosome(2),
                (i * 1000) as u32,
                i as f64 * 0.01, // 1 cM = 0.01 Morgans
                'A',
                'G',
            ));
        }

        let blocks = BlockDefinition::default().assign_blocks(&snps);

        // Should have 2 blocks (one per chromosome, each <5 cM)
        assert_eq!(blocks.n_blocks, 2);
        assert_eq!(blocks.block_of(0), 0); // First chr
        assert_eq!(blocks.block_of(4), 0);
        assert_eq!(blocks.block_of(5), 1); // Second chr
        assert_eq!(blocks.block_of(9), 1);
    }

    #[test]
    fn test_empty_snps() {
        let snps: Vec<Snp> = Vec::new();
        let blocks = BlockDefinition::default().assign_blocks(&snps);

        assert_eq!(blocks.n_blocks, 0);
        assert!(blocks.snp_to_block.is_empty());
    }

    #[test]
    fn test_custom_block_size() {
        // 10 SNPs at 1 cM spacing with 2 cM blocks
        let snps = make_snps_same_chr(10, 1.0);
        let blocks = BlockDefinition::new(0.02).assign_blocks(&snps); // 2 cM = 0.02 M blocks

        // Allow tolerance for floating point precision at boundaries
        assert!(
            blocks.n_blocks >= 4 && blocks.n_blocks <= 6,
            "Expected ~5 blocks, got {}",
            blocks.n_blocks
        );
    }

    #[test]
    fn test_block_range() {
        let snps = make_snps_same_chr(100, 1.0);
        let blocks = BlockDefinition::default().assign_blocks(&snps);

        // First block should be 0..5
        assert_eq!(blocks.block_range(0), 0..5);
        // Second block should be 5..10
        assert_eq!(blocks.block_range(1), 5..10);
    }
}
