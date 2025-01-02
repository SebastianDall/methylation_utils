use ahash::AHashMap;
use anyhow::{anyhow, Context};
use log::debug;
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

        for record_result in self.reader.records() {
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
                debug!("Current contig id in line: {}", &contig_id);

                debug!(
                    "Current contig being added: {}",
                    self.current_contig
                        .as_ref()
                        .map(|c| c.id.to_string())
                        .unwrap_or("None".to_string())
                );
                if let Some(old_contig) = self.current_contig.take() {
                    debug!("Adding contig to builder");
                    if let Err(e) = builder.add_contig(old_contig) {
                        return Some(Err(e));
                    }

                    debug!(
                        "Contigs loaded in batch: {}. Batch size is: {}",
                        self.contigs_loaded_in_batch, self.batch_size
                    );
                    if self.contigs_loaded_in_batch == self.batch_size {
                        self.contigs_loaded_in_batch = 0;
                        return Some(Ok(builder.build()));
                    }
                };

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

#[cfg(test)]
mod tests {
    use crate::data::methylation::MethylationCoverage;

    use super::*;
    use std::{
        fs::File,
        io::{BufReader, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_batch_loading() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let mut assembly = AHashMap::new();
        assembly.insert(
            "contig_3".to_string(),
            Contig::new("contig_3".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = BatchLoader::new(reader, assembly, 1, 1);

        for ws in batch_loader {
            let workspace = ws?.get_workspace();
            assert_eq!(workspace.len(), 1);

            assert_eq!(
                workspace
                    .get("contig_3")
                    .unwrap()
                    .get_methylated_positions(
                        &[6],
                        methylome::Strand::Positive,
                        methylome::ModType::SixMA
                    )
                    .into_iter()
                    .map(|rec| rec.unwrap().clone())
                    .collect::<Vec<MethylationCoverage>>(),
                vec![MethylationCoverage::new(15, 15).unwrap()]
            );
        }

        Ok(())
    }

    #[test]
    fn test_multi_batches() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let mut assembly = AHashMap::new();
        assembly.insert(
            "contig_3".to_string(),
            Contig::new("contig_3".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        assembly.insert(
            "contig_4".to_string(),
            Contig::new("contig_4".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = BatchLoader::new(reader, assembly, 1, 1);

        let mut num_batches = 0;
        for _ws in batch_loader {
            num_batches += 1;
        }

        assert_eq!(num_batches, 2);

        Ok(())
    }

    #[test]
    fn test_contig_missing_error() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let assembly = AHashMap::new();
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = BatchLoader::new(reader, assembly, 2, 1);

        for ws in batch_loader {
            assert!(ws.is_err());
        }

        Ok(())
    }
}
