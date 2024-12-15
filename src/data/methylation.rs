use anyhow::{bail, Result};

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct MethylationCoverage {
    n_modified: u32,
    n_valid_cov: u32,
}

impl MethylationCoverage {
    pub fn new(n_modified: u32, n_valid_cov: u32) -> Result<Self> {
        if n_modified > n_valid_cov {
            bail!(
                "Invalid coverage: n_valid_cov ({}) cannot be less than n_modified ({})",
                n_valid_cov,
                n_modified
            )
        }

        Ok(Self {
            n_modified,
            n_valid_cov,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
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
