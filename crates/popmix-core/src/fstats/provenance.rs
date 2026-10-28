//! Where f2 blocks came from (TASKS.md H3): a JSON object the code that
//! makes them fills in (`popmix_io::extract`), saved in an f2 directory's
//! `metadata.json` and kept through population subsets and block
//! resampling. Blocks made otherwise (from frequencies or block means)
//! have none.

use super::f2_blocks::F2Blocks;

impl F2Blocks {
    /// How the blocks were made, if known.
    pub fn provenance(&self) -> Option<&serde_json::Value> {
        self.provenance.as_ref()
    }

    /// The same blocks, recorded as made as `provenance` says.
    #[must_use]
    pub fn with_provenance(mut self, provenance: Option<serde_json::Value>) -> Self {
        self.provenance = provenance;
        self
    }
}
