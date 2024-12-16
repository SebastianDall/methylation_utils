pub mod contig;
pub mod methylation;

use anyhow::{bail, Context, Result};
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use log::{error, info};
use methylation::MethylationCoverage;
use methylome::{ModType, Strand};

use crate::data::contig::Contig;
use std::{
    collections::HashMap,
    fmt::Write,
    fs::File,
    io::{BufRead, BufReader},
};

pub struct GenomeWorkspace {
    pub contigs: HashMap<String, Contig>,
}

impl GenomeWorkspace {
    pub fn new() -> Self {
        Self {
            contigs: HashMap::new(),
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
        let file = File::open(&pileup_path)
            .with_context(|| format!("Error opening file: {:?}", &pileup_path.as_ref()))?;

        let reader = BufReader::new(&file);
        let total_lines = reader.lines().count();

        let pb = ProgressBar::new(total_lines as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len}",
            )
            .unwrap()
            .with_key("eta", |state: &ProgressState, w: &mut dyn Write| {
                write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap()
            })
            .progress_chars("#>-"),
        );

        let file = File::open(&pileup_path)
            .with_context(|| format!("Error re-opening file: {:?}", &pileup_path.as_ref()))?;

        let reader = BufReader::new(&file);

        for (line_num, line_result) in reader.lines().enumerate() {
            pb.inc(1);

            let line = line_result.with_context(|| {
                format!(
                    "Could not read line '{}' in file: {:?}",
                    line_num + 1,
                    &pileup_path.as_ref()
                )
            })?;

            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() != 18 {
                pb.finish_with_message("Pileup processing complete with errors");
                bail!(
                    "Malformed line on line '{}' in file: {:?}",
                    line_num + 1,
                    &pileup_path.as_ref()
                )
            }
            let n_modified: u32 = fields[11].parse().with_context(|| {
                format!(
                    "Invalid N_modified value: {}. Occured at line {}",
                    fields[11],
                    line_num + 1
                )
            })?;
            let n_valid_cov: u32 = fields[9].parse().with_context(|| {
                format!(
                    "Invalid N_valid_cov value: {}. Occured at line {}",
                    fields[9],
                    line_num + 1
                )
            })?;

            // Skip entries with zero coverage
            if n_valid_cov < min_valid_read_coverage {
                continue;
            }

            let methylation = MethylationCoverage::new(n_modified, n_valid_cov)?;

            let contig_id = fields[0];
            let position: usize = fields[1].parse().with_context(|| {
                format!("Could not parse position to int on line: {}", line_num + 1)
            })?;
            let strand = Strand::from_str(fields[5])?;
            let mod_type = ModType::from_str(fields[3])?;

            if let Some(contig_entry) = self.get_mut_contig(&contig_id) {
                contig_entry.add_methylation(position, strand, mod_type, methylation)?;
            } else {
                error!(
                    "Warning: Contig: '{}' found in pileup, but not in assembly",
                    contig_id
                );
            }
        }
        pb.finish_with_message("Finished loading methylation from pileup.");
        info!("Pileup processing succesfull");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use std::fs::File;
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