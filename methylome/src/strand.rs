use std::fmt::Display;

use anyhow::{bail, Result};

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
    pub fn from_str(strand: &str) -> Result<Self> {
        match strand {
            "+" => Ok(Strand::Positive),
            "-" => Ok(Strand::Negative),
            _ => bail!("Could not parse '{}' to Strand", strand),
        }
    }

    pub fn to_string(&self) -> String {
        match self {
            Strand::Positive => "+".to_string(),
            Strand::Negative => "-".to_string(),
        }
    }
}
