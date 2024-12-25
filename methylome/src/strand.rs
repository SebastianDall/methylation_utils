use std::{fmt::Display, str::FromStr};

use anyhow::{bail, Result};

/// Represents the DNA strand of reference.
#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Positive,
    Negative,
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

impl Strand {
    pub fn to_string(&self) -> String {
        match self {
            Strand::Positive => "+".to_string(),
            Strand::Negative => "-".to_string(),
        }
    }
}

/// Parses a &str to the Strand enum type.
/// Should be either:
/// - +: Positive
/// - -: Negative
///
/// # Examples
/// ```
/// use methylome::Strand;
///
/// let strand = "+".parse::<Strand>().unwrap();
/// assert_eq!(strand, Strand::Positive);
///
/// let invalid_strand = "p".parse::<Strand>();
/// assert!(invalid_strand.is_err());
/// ```
impl FromStr for Strand {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Positive),
            "-" => Ok(Strand::Negative),
            _ => bail!("Could not parse '{}' to Strand", s),
        }
    }
}
