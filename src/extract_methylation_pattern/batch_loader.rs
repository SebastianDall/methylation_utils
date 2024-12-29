use ahash::HashMap;
use std::io::{BufRead, BufReader};

use crate::data::{contig::Contig, GenomeWorkspace, GenomeWorkspaceBuilder};

pub struct BatchLoader<R> {
    reader: csv::Reader<R>,
    assembly: HashMap<String, Contig>,
    batch_size: usize,

    current_contig_id: Option<String>,
    current_contig: Option<Contig>,
    contigs_loaded_in_batch: usize,
}

impl<R: BufRead> BatchLoader<R> {
    pub fn new(reader: R, assembly: HashMap<String, Contig>, batch_size: usize) -> Self {
        let rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .flexible(false)
            .from_reader(reader);

        BatchLoader {
            reader: rdr,
            assembly: assembly,
            batch_size: batch_size,
            current_contig_id: None,
            current_contig: None,
            contigs_loaded_in_batch: 0,
        }
    }
}

impl<R: BufRead> Iterator for BatchLoader<R> {
    type Item = GenomeWorkspace;

    fn next(&mut self) -> Option<GenomeWorkspace> {
        let mut builder = GenomeWorkspaceBuilder::new();

        let n_valid_cov: u32 = record
            .get(9)
            .ok_or_else(|| anyhow!("Missing n_valid_coverage field"))?
            .parse()
            .map_err(|_| anyhow!("Invalid coverage number."))?;
        if n_valid_cov < args.min_valid_read_coverage {
            continue;
        }

        let contig_id = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing contig field"))?
            .to_string();
    }
}
