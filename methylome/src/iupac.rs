use std::fmt::Display;

use anyhow::bail;

// INFO: This is taken from modkit https://github.com/nanoporetech/modkit?tab=License-1-ov-file
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IupacBase {
    A,
    T,
    G,
    C,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    N,
}

impl Display for IupacBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

impl IupacBase {
    pub fn parse_char(base: char) -> anyhow::Result<Self> {
        let iupac_base = match base {
            'A' => Ok(Self::A),
            'T' => Ok(Self::T),
            'G' => Ok(Self::G),
            'C' => Ok(Self::C),
            'R' => Ok(Self::R),
            'Y' => Ok(Self::Y),
            'S' => Ok(Self::S),
            'W' => Ok(Self::W),
            'K' => Ok(Self::K),
            'M' => Ok(Self::M),
            'B' => Ok(Self::B),
            'D' => Ok(Self::D),
            'H' => Ok(Self::H),
            'V' => Ok(Self::V),
            'N' => Ok(Self::N),
            _ => bail!("Not a defined Iupac base: {base}"),
        };
        iupac_base
    }

    pub fn to_string(&self) -> String {
        match self {
            IupacBase::A => "A".to_string(),
            IupacBase::T => "T".to_string(),
            IupacBase::G => "G".to_string(),
            IupacBase::C => "C".to_string(),
            IupacBase::R => "R".to_string(),
            IupacBase::Y => "Y".to_string(),
            IupacBase::S => "S".to_string(),
            IupacBase::W => "W".to_string(),
            IupacBase::K => "K".to_string(),
            IupacBase::M => "M".to_string(),
            IupacBase::B => "B".to_string(),
            IupacBase::D => "D".to_string(),
            IupacBase::H => "H".to_string(),
            IupacBase::V => "V".to_string(),
            IupacBase::N => "N".to_string(),
        }
    }

    pub fn to_complement_base(base: &IupacBase) -> Self {
        match base {
            IupacBase::A => IupacBase::T,
            IupacBase::T => IupacBase::A,
            IupacBase::G => IupacBase::C,
            IupacBase::C => IupacBase::G,
            IupacBase::R => IupacBase::Y,
            IupacBase::Y => IupacBase::R,
            IupacBase::S => IupacBase::S,
            IupacBase::W => IupacBase::W,
            IupacBase::K => IupacBase::M,
            IupacBase::M => IupacBase::K,
            IupacBase::B => IupacBase::V,
            IupacBase::D => IupacBase::H,
            IupacBase::H => IupacBase::D,
            IupacBase::V => IupacBase::B,
            IupacBase::N => IupacBase::N,
        }
    }

    pub fn to_regex(&self) -> &str {
        match self {
            IupacBase::A => "A",
            IupacBase::T => "T",
            IupacBase::G => "G",
            IupacBase::C => "C",
            IupacBase::R => "[AG]",
            IupacBase::Y => "[CT]",
            IupacBase::S => "[CG]",
            IupacBase::W => "[AT]",
            IupacBase::K => "[GT]",
            IupacBase::M => "[AC]",
            IupacBase::B => "[CGT]",
            IupacBase::D => "[AGT]",
            IupacBase::H => "[ACT]",
            IupacBase::V => "[ACG]",
            IupacBase::N => ".",
        }
    }
}
