use anyhow::{bail, Result};
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModType {
    SixMA,
    FiveMC,
    FourMC,
}

impl ModType {
    pub fn from_str(mod_type: &str) -> Result<Self> {
        match mod_type {
            "a" => Ok(ModType::SixMA),
            "m" => Ok(ModType::FiveMC),
            "21839" => Ok(ModType::FourMC),
            _ => bail!("Unsupported mod type: {}", mod_type),
        }
    }

    pub fn to_pileup_code(&self) -> &'static str {
        match self {
            ModType::SixMA => "a",
            ModType::FiveMC => "m",
            ModType::FourMC => "21839",
        }
    }
}

impl fmt::Display for ModType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ModType::SixMA => write!(f, "6mA (a)"),
            ModType::FiveMC => write!(f, "5mC (m)"),
            ModType::FourMC => write!(f, "4mC (21839)"),
        }
    }
}
