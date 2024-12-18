pub mod contig;
pub mod methylation;

use crate::data::contig::Contig;
use ahash::AHashMap;
use anyhow::{anyhow, bail, Context, Result};
use bytesize::ByteSize;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use log::{error, info};
use methylation::MethylationCoverage;
use methylome::{ModType, Strand};
use std::{fmt::Write, fs::read, io::Cursor};

struct MethylationRecord {
    contig: String,
    position: usize,
    strand: Strand,
    mod_type: ModType,
    methylation: MethylationCoverage,
}

pub struct GenomeWorkspace {
    pub contigs: AHashMap<String, Contig>,
}

impl GenomeWorkspace {
    pub fn new() -> Self {
        Self {
            contigs: AHashMap::new(),
        }
    }

    pub fn add_contig(&mut self, contig: Contig) -> Result<()> {
        if self.contigs.contains_key(&contig.id) {
            bail!("Key error: '{}' already inserted", &contig.id)
        }

        self.contigs.insert(contig.id.clone(), contig);
        Ok(())
    }

    pub fn prune_empty_contigs(&mut self) {
        self.contigs
            .retain(|_id, contig| !contig.methylated_positions.is_empty());
    }

    pub fn get_mut_contig(&mut self, id: &str) -> Option<&mut Contig> {
        self.contigs.get_mut(id)
    }

    pub fn populate_methylation_from_pileup<P: AsRef<std::path::Path>>(
        &mut self,
        pileup_path: P,
        min_valid_read_coverage: u32,
    ) -> Result<()> {
        let pileup_path_ref = pileup_path.as_ref();

        let data = read(pileup_path_ref)
            .with_context(|| format!("Error reading entire file: {:?}", pileup_path_ref))?;

        let file_size = data.len() as u64;
        let human_readable_size = ByteSize::b(file_size).to_string();

        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes:>8}/{total_bytes:>8} ({percent}%)",
            )
            .unwrap()
            .with_key("bytes", |state: &ProgressState, w: &mut dyn Write| {
                // Convert processed bytes to human-readable format
                write!(w, "{}", ByteSize::b(state.pos())).unwrap()
            })
                .with_key("total_bytes", {
                let human_readable_size = human_readable_size.clone();
                move |_state: &ProgressState, w: &mut dyn std::fmt::Write| {
                    write!(w, "{}", human_readable_size).unwrap()
                }
            })
            .progress_chars("#>-"),
        );

        let cursor = Cursor::new(data);
        let mut rdr = ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(cursor);
        let mut record = csv::StringRecord::new();

        let mut methylation_records: Vec<MethylationRecord> = Vec::new();

        let mut lines_processed = 0;
        // for (line_num, result) in rdr.records().enumerate() {
        while rdr.read_record(&mut record)? {
            lines_processed += 1;

            if lines_processed % 100_000 == 0 {
                let bytes_read = rdr.position().byte();
                pb.set_position(bytes_read as u64);
            }

            let n_valid_cov_str = record
                .get(9)
                .with_context(|| anyhow!("Missing contig at line {}", lines_processed))?;
            let n_valid_cov = n_valid_cov_str.parse().with_context(|| {
                format!(
                    "Invalid N_valid_cov value: {}. Occured at line {}",
                    n_valid_cov_str, lines_processed
                )
            })?;
            // Skip entries with zero coverage
            if n_valid_cov < min_valid_read_coverage {
                continue;
            }

            let contig_id = record
                .get(0)
                .ok_or_else(|| anyhow!("Missing contig at line {}", lines_processed))?
                .to_string();

            let position_str = record
                .get(1)
                .ok_or_else(|| anyhow!("Missing position at line {}", lines_processed))?;
            let position: usize = position_str
                .parse()
                .with_context(|| format!("Invalid position at line {}", lines_processed))?;

            let mod_type_str = record
                .get(3)
                .ok_or_else(|| anyhow!("Missing mod_type at line {}", lines_processed))?;
            let mod_type = ModType::from_str(mod_type_str)?;

            let strand_str = record
                .get(5)
                .ok_or_else(|| anyhow!("Missing strand at line {}", lines_processed))?;
            let strand = Strand::from_str(strand_str)?;

            let n_modified_str = record
                .get(11)
                .ok_or_else(|| anyhow!("Missing n_modified at line {}", lines_processed))?;
            let n_modified: u32 = n_modified_str
                .parse()
                .with_context(|| format!("Invalid n_modified at line {}", lines_processed))?;

            let methylation = MethylationCoverage::new(n_modified, n_valid_cov)?;

            methylation_records.push(MethylationRecord {
                contig: contig_id,
                position,
                strand,
                mod_type,
                methylation,
            });
        }
        pb.finish_with_message("Finished reading pileup into memory.\n");

        for rec in methylation_records {
            if let Some(contig_entry) = self.get_mut_contig(&rec.contig) {
                contig_entry.add_methylation(
                    rec.position,
                    rec.strand,
                    rec.mod_type,
                    rec.methylation,
                )?;
            } else {
                error!(
                    "Warning: Contig: '{}' found in pileup, but not in assembly",
                    rec.contig
                );
            }
        }

        info!("Pileup processing succesfull");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_strand_from_str() -> Result<()> {
        // Mock pileup data lines
        let pileup_data = vec![
            "contig_3\t0\t1\tm\t133\t-\t0\t1\t255,0,0\t133\t0.00\t0\t133\t0\t0\t6\t0\t0",
            "contig_3\t1\t2\ta\t174\t+\t1\t2\t255,0,0\t174\t1.72\t3\t171\t0\t0\t3\t0\t0",
        ];

        // Expected results for the strand column
        let expected_strands = vec![Strand::Negative, Strand::Positive];

        // Iterate through pileup data and validate strand parsing
        for (line, &expected_strand) in pileup_data.iter().zip(expected_strands.iter()) {
            let fields: Vec<&str> = line.split('\t').collect();
            let strand_field = fields[5]; // Extract strand field
            let strand = Strand::from_str(strand_field)?;

            // Assert that the parsed strand matches the expected value
            assert_eq!(strand, expected_strand);
        }

        Ok(())
    }

    #[test]
    fn test_populate_methylation() -> Result<()> {
        let mut workspace = GenomeWorkspace::new();

        // Add a mock contig to the workspace
        workspace.add_contig(Contig::new("contig_3".to_string(), "ATCG".to_string()))?;

        // Create a temporary pileup file
        let mut pileup_file = NamedTempFile::new()?;
        writeln!(
            pileup_file,
            "contig_3\t0\t1\tm\t133\t-\t0\t1\t255,0,0\t133\t0.00\t10\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t1\t2\ta\t174\t+\t1\t2\t255,0,0\t174\t1.72\t5\t169\t0\t0\t3\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t2\t3\ta\t172\t+\t2\t3\t255,0,0\t0\t0.00\t0\t0\t0\t0\t0\t0\t0" // Zero coverage, should be skipped
        )?;

        // Populate methylation data
        workspace.populate_methylation_from_pileup(pileup_file.path(), 1)?;

        // Get the contig
        let contig = workspace.get_mut_contig("contig_3").unwrap();

        // Check that methylated_positions are correctly populated
        assert_eq!(contig.methylated_positions.len(), 2);

        // Validate individual positions
        let pos_0 = contig
            .methylated_positions
            .get(&(0, Strand::Negative, ModType::FiveMC));

        let expected_methylation = MethylationCoverage::new(10, 133).unwrap();

        assert!(pos_0.is_some());
        assert_eq!(pos_0, Some(&expected_methylation));

        let pos_1 = contig
            .methylated_positions
            .get(&(1, Strand::Positive, ModType::SixMA));
        let expected_methylation = MethylationCoverage::new(5, 174).unwrap();
        assert!(pos_1.is_some());
        assert_eq!(pos_1, Some(&expected_methylation));

        // Ensure position with zero coverage is skipped
        let pos_2 = contig
            .methylated_positions
            .get(&(2, Strand::Positive, ModType::SixMA));
        assert!(pos_2.is_none());

        Ok(())
    }

    #[test]
    fn test_populate_methylation_missing_contig() -> Result<()> {
        let mut workspace = GenomeWorkspace::new();

        // Add a mock contig that doesn't match the pileup
        workspace.add_contig(Contig::new("contig_1".to_string(), "ATCG".to_string()))?;

        // Create a temporary pileup file
        let mut pileup_file = NamedTempFile::new()?;
        writeln!(
            pileup_file,
            "contig_3\t0\t1\tm\t133\t-\t0\t1\t255,0,0\t133\t0.00\t10\t123\t0\t0\t6\t0\t0"
        )?;

        // Populate methylation data
        workspace.populate_methylation_from_pileup(pileup_file.path(), 1)?;

        // Ensure the workspace is unchanged since the contig doesn't match
        let contig = workspace.get_mut_contig("contig_1").unwrap();
        assert!(contig.methylated_positions.is_empty());

        Ok(())
    }
}
