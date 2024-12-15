use std::{collections::HashMap, fmt::Display};

use anyhow::{bail, Result};

mod data;

use crate::ModType;
use data::methylation::*;

pub struct Contig {
    pub id: String,
    pub sequence: String,
    sequence_len: usize,
    methylated_positions: HashMap<(usize, Strand, ModType), MethylationCoverage>,
}

impl Contig {
    pub fn new(id: String, sequence: String) -> Self {
        let sequence_length = sequence.len();

        Self {
            id,
            sequence,
            sequence_len: sequence_length,
            methylated_positions: HashMap::new(),
        }
    }

    pub fn add_methylation(
        &mut self,
        position: usize,
        strand: Strand,
        mod_type: ModType,
        meth_coverage: MethylationCoverage,
    ) -> Result<()> {
        if position as usize >= self.sequence_len {
            bail!("Position out of bounds for '{}': Cannot insert key position ({}) longer than contig length ({})!", self.id, position, self.sequence_len)
        }

        let key = (position, strand.clone(), mod_type.clone());

        if self.methylated_positions.contains_key(&key) {
            bail!("Methylation record already store for: {} - strand ({}) - modification type ({}) - position '{}'",self.id, strand,mod_type, position)
        }

        self.methylated_positions.insert(key, meth_coverage);
        Ok(())
    }

    pub fn get_methylated_positions(
        &self,
        positions: &[usize],
        strand: Strand,
        mod_type: ModType,
    ) -> Vec<Option<&MethylationCoverage>> {
        positions
            .iter()
            .map(|&pos| self.methylated_positions.get(&(pos, strand, mod_type)))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_construction() {
        let mut contig = Contig::new("contig_1".to_string(), "TGGACGATCCCGATC".to_string());

        let meth_record1 = MethylationCoverage::new(1, 1).unwrap();
        let meth_record2 = MethylationCoverage::new(2, 2).unwrap();
        let meth_record3 = MethylationCoverage::new(3, 3).unwrap();

        // Insert 6mA records
        contig
            .add_methylation(6, Strand::Positive, ModType::SixMA, meth_record1.clone())
            .unwrap();
        contig
            .add_methylation(12, Strand::Positive, ModType::SixMA, meth_record1.clone())
            .unwrap();
        contig
            .add_methylation(13, Strand::Negative, ModType::SixMA, meth_record1.clone())
            .unwrap();

        // Insert 5mC record
        contig
            .add_methylation(8, Strand::Positive, ModType::FiveMC, meth_record3)
            .unwrap();

        // Insert unused record that should not be returned
        contig
            .add_methylation(6, Strand::Positive, ModType::FiveMC, meth_record2.clone())
            .unwrap();

        let positions: Vec<usize> = vec![6, 12];

        let meth_records =
            contig.get_methylated_positions(&positions, Strand::Positive, ModType::SixMA);

        // Ensure records match the expected values
        let expected = vec![
            Some(&MethylationCoverage {
                n_modified: 1,
                n_valid_cov: 1,
            }),
            Some(&MethylationCoverage {
                n_modified: 1,
                n_valid_cov: 1,
            }),
        ];

        assert_eq!(meth_records, expected);

        let meth_records = contig.get_methylated_positions(&[13], Strand::Negative, ModType::SixMA);
        let expected = vec![Some(&MethylationCoverage {
            n_modified: 1,
            n_valid_cov: 1,
        })];

        assert_eq!(meth_records, expected);

        let meth_records = contig.get_methylated_positions(&[8], Strand::Positive, ModType::FiveMC);
        assert_eq!(
            meth_records,
            vec![Some(&MethylationCoverage {
                n_modified: 3,
                n_valid_cov: 3
            })]
        )
    }

    #[test]
    fn test_out_of_bounds_record() {
        let mut contig = Contig::new("1".to_string(), "GATC".to_string());

        let result = contig.add_methylation(
            3,
            Strand::Positive,
            ModType::SixMA,
            MethylationCoverage {
                n_modified: 1,
                n_valid_cov: 1,
            },
        );

        assert!(result.is_ok());
    }

    #[test]
    fn test_methylation_coverage_valid() -> Result<()> {
        // Test valid inputs
        let coverage = MethylationCoverage::new(5, 10)?;
        assert_eq!(coverage.n_modified, 5);
        assert_eq!(coverage.n_valid_cov, 10);

        let coverage = MethylationCoverage::new(0, 0)?;
        assert_eq!(coverage.n_modified, 0);
        assert_eq!(coverage.n_valid_cov, 0);

        Ok(())
    }

    #[test]
    fn test_methylation_coverage_invalid() {
        // Test invalid input: n_valid_cov < n_modified
        let result = MethylationCoverage::new(10, 5);

        assert!(result.is_err());
        if let Err(e) = result {
            assert_eq!(
                e.to_string(),
                "Invalid coverage: n_valid_cov (5) cannot be less than n_modified (10)"
            );
        }
    }
}
