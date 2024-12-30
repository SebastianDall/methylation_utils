use ahash::AHashMap;
use anyhow::{anyhow, Context};
use std::io::BufRead;

use crate::data::{contig::Contig, GenomeWorkspace, GenomeWorkspaceBuilder};

use super::parse_to_methylation_record;

pub struct BatchLoader<R> {
    reader: csv::Reader<R>,
    assembly: AHashMap<String, Contig>,
    batch_size: usize,
    min_valid_read_coverage: u32,

    current_contig_id: Option<String>,
    current_contig: Option<Contig>,
    contigs_loaded_in_batch: usize,
}

impl<R: BufRead> BatchLoader<R> {
    pub fn new(
        reader: R,
        assembly: AHashMap<String, Contig>,
        batch_size: usize,
        min_valid_read_coverage: u32,
    ) -> Self {
        let rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .flexible(false)
            .from_reader(reader);

        let size = if batch_size == 0 { 1 } else { batch_size };

        BatchLoader {
            reader: rdr,
            assembly,
            batch_size: size,
            min_valid_read_coverage,
            current_contig_id: None,
            current_contig: None,
            contigs_loaded_in_batch: 0,
        }
    }
}

impl<R: BufRead> Iterator for BatchLoader<R> {
    type Item = Result<GenomeWorkspace, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut builder = GenomeWorkspaceBuilder::new();

        while let Some(record_result) = self.reader.records().next() {
            let record = match record_result.context("Failed to read pileup record") {
                Ok(r) => r,
                Err(e) => return Some(Err(e)),
            };

            let n_valid_cov: u32 = match record.get(9) {
                Some(f) => match f.parse() {
                    Ok(v) => v,
                    Err(_) => return Some(Err(anyhow!("Invalid coverage number"))),
                },
                None => return Some(Err(anyhow!("Invalid coverage field"))),
            };

            if n_valid_cov < self.min_valid_read_coverage {
                continue;
            }

            let contig_id = match record.get(0) {
                Some(id) => id.to_string(),
                None => return Some(Err(anyhow!("Missing contig id field"))),
            };

            if Some(&contig_id) != self.current_contig_id.as_ref() {
                if let Some(old_contig) = self.current_contig.take() {
                    if let Err(e) = builder.add_contig(old_contig) {
                        return Some(Err(e));
                    }

                    if self.contigs_loaded_in_batch > self.batch_size {
                        return Some(Ok(builder.build()));
                    }

                    self.current_contig_id = Some(contig_id.clone());
                    self.contigs_loaded_in_batch += 1;

                    let contig_key = self.current_contig_id.as_ref().unwrap();
                    let new_contig = match self.assembly.get(contig_key) {
                        Some(c) => c.clone(),
                        None => {
                            return Some(Err(anyhow!(
                                "Contig '{}' not found in assembly",
                                contig_key
                            )))
                        }
                    };
                    self.current_contig = Some(new_contig);
                }
                let meth = match parse_to_methylation_record(contig_id, n_valid_cov, &record) {
                    Ok(m) => m,
                    Err(e) => return Some(Err(e)),
                };
                if let Some(ref mut c) = self.current_contig {
                    if let Err(e) = c.add_methylation_record(meth) {
                        return Some(Err(e));
                    }
                }
            }
        }
        if let Some(last) = self.current_contig.take() {
            builder.add_contig(last).ok()?;
        }

        let workspace = builder.build();
        if workspace.is_empty() {
            None
        } else {
            Some(Ok(workspace))
        }
    }
}
