use anyhow::bail;
use phf::phf_map;
use regex::Regex;
use std::fmt;

static IUPAC_MAPPING: phf::Map<&'static str, &'static str> = phf_map! {
    "R" => "AG",
    "Y" => "CT",
    "S" => "CG",
    "W" => "AT",
    "K" => "GT",
    "M" => "AC",
    "B" => "CGT",
    "D" => "AGT",
    "H" => "ACT",
    "V" => "ACG",
    "N" => ".",
};

static COMPLEMENT_BASE: phf::Map<&'static str, &'static str> = phf_map! {
    "A" => "T",
    "T" => "A",
    "G" => "C",
    "C" => "G",
    "R" => "Y",
    "Y" => "R",
    "S" => "S",
    "W" => "W",
    "K" => "M",
    "M" => "K",
    "B" => "V",
    "D" => "H",
    "H" => "D",
    "V" => "B",
    "." => ".",
    "[" => "]",
    "]" => "[",
    "N" => "N",
};

// INFO: This is taken from modkit https://github.com/nanoporetech/modkit?tab=License-1-ov-file
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

impl IupacBase {
    pub fn parse_char(base: char) -> anyhow::Result<Self> {
        match base {
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
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModType {
    SixMA,
    FiveMC,
    FourMC,
}

impl ModType {
    pub fn from_str(mod_type: &str) -> Result<Self, String> {
        match mod_type {
            "a" => Ok(ModType::SixMA),
            "m" => Ok(ModType::FiveMC),
            "21839" => Ok(ModType::FourMC),
            _ => Err(format!("Unsupported mod type: {}", mod_type)),
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Motif {
    pub sequence: String,
    pub mod_type: ModType,
    pub mod_position: u8,
}

impl Motif {
    pub fn new(sequence: &str, mod_type: &str, mod_position: u8) -> Result<Self, String> {
        let mod_type = ModType::from_str(mod_type)?;

        for b in sequence.chars() {
            IupacBase::parse_char(b).map_err(|_| {
                format!(
                    "Base '{}' in sequence '{}' is not a valid IUPAC code",
                    b, sequence
                )
            })?;
        }

        if mod_position as usize > sequence.len() - 1 {
            return Err(format!(
                "mod_position {} is out of bounds for sequence of length {}. Note mod_position is 0-indexed.",
                mod_position,
                sequence.len()
            ));
        }

        let base_at_position = sequence.chars().nth(mod_position as usize).unwrap(); // Safe as we have checked bounds above
        match mod_type {
            ModType::SixMA => {
                if base_at_position != 'A' {
                    return Err(format!(
                        "mod_position {} points to base '{}' which is invalid for 6mA",
                        mod_position, base_at_position
                    ));
                }
            }
            ModType::FiveMC | ModType::FourMC => {
                if base_at_position != 'C' {
                    return Err(format!(
                        "mod_position {} points to base '{}' which is invalid for {} modification type.",
                        mod_position, base_at_position, mod_type
                    ));
                }
            }
        }

        Ok(Self {
            sequence: sequence.to_string(),
            mod_type,
            mod_position,
        })
    }

    pub fn reverse_complement(&self) -> Self {
        Self {
            // sequence: (&self.sequence.chars().rev().collect::<String>()).to_string(),
            sequence: self
                .sequence
                .chars()
                .rev()
                .map(|c| match COMPLEMENT_BASE.get(&c.to_string().as_str()) {
                    Some(s) => s.to_string(),
                    None => c.to_string(),
                })
                .collect(),
            mod_type: self.mod_type.clone(),
            mod_position: self.sequence.len() as u8 - self.mod_position - 1,
        }
    }

    pub fn to_regex(&self) -> String {
        self.sequence
            .chars()
            .map(|c| match IUPAC_MAPPING.get(&c.to_string().as_str()) {
                Some(s) => {
                    if *s == "." {
                        ".".to_string()
                    } else {
                        format!("[{}]", s)
                    }
                }
                None => c.to_string(),
            })
            .collect()
    }
}

pub fn find_motif_indices_in_contig(contig: &str, motif: &Motif) -> Vec<u32> {
    let regex_str = motif.to_regex();
    let re = Regex::new(&regex_str).expect("Expected regex pattern");

    let indices = re
        .find_iter(contig)
        .map(|m| m.start() as u32 + motif.mod_position as u32)
        .collect();

    indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_creation() {
        let motif = Motif::new("GATC", "a", 1).unwrap();
        assert_eq!(motif.sequence, "GATC");
        assert_eq!(motif.mod_type, ModType::SixMA);
        assert_eq!(motif.mod_position, 1);
    }

    #[test]
    fn test_out_of_bounds() {
        let result = Motif::new("GATC", "m", 4);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "mod_position 4 is out of bounds for sequence of length 4. Note mod_position is 0-indexed."
        );
    }

    #[test]
    fn test_unidentified_motif_type() {
        let result = Motif::new("GATC", "d", 1);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Unsupported mod type: d");
    }

    #[test]
    fn test_invalid_mod_position_base() {
        let result = Motif::new("ATCG", "m", 3); // 'G' is invalid for 5mC
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "mod_position 3 points to base 'G' which is invalid for 5mC (m) modification type."
        );
    }

    #[test]
    fn test_invalid_iupac_base() {
        let result = Motif::new("ATZG", "a", 0); // 'G' is invalid for 5mC
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Base 'Z' in sequence 'ATZG' is not a valid IUPAC code"
        );
    }

    #[test]
    fn test_motif_reverse_complement() {
        let motif1 = Motif::new("GATC", "m", 3).unwrap();
        let motif2 = Motif::new("TCCCG", "m", 1).unwrap();
        let motif3 = Motif::new("RGATCY", "a", 2).unwrap();
        assert_eq!(motif1.reverse_complement().sequence, "GATC");
        assert_eq!(motif2.reverse_complement().sequence, "CGGGA");
        assert_eq!(motif3.reverse_complement().sequence, "RGATCY");
        assert_eq!(
            motif1.reverse_complement().mod_type,
            ModType::from_str("m").unwrap()
        );
        assert_eq!(
            motif2.reverse_complement().mod_type,
            ModType::from_str("m").unwrap()
        );
        assert_eq!(
            motif3.reverse_complement().mod_type,
            ModType::from_str("a").unwrap()
        );
        assert_eq!(motif1.reverse_complement().mod_position, 0);
        assert_eq!(motif2.reverse_complement().mod_position, 3);
        assert_eq!(motif3.reverse_complement().mod_position, 3);
    }

    #[test]
    fn test_to_regex() {
        let motif1 = Motif::new("GATC", "m", 3).unwrap();
        let motif2 = Motif::new("RGATCY", "m", 4).unwrap();

        assert_eq!(motif1.to_regex(), "GATC");
        assert_eq!(motif2.to_regex(), "[AG]GATC[CT]");
    }

    #[test]
    fn test_find_motif_indices_in_contig() {
        let contig = "GGATCTCCATGATC".to_string();
        let contig2 = "TGGACGATCCCGATC".to_string();
        let motif1 = Motif::new("GATC", "m", 3).unwrap();
        let motif2 = Motif::new("RGATCY", "m", 4).unwrap();
        let motif3 = Motif::new("GATC", "a", 1).unwrap();
        let motif4 = Motif::new("GGANNNTCC", "a", 2).unwrap();

        println!("{}", &motif4.to_regex());
        assert_eq!(find_motif_indices_in_contig(&contig, &motif1), vec![4, 13]);
        assert_eq!(find_motif_indices_in_contig(&contig, &motif2), vec![4]);

        assert_eq!(find_motif_indices_in_contig(&contig2, &motif3), vec![6, 12]);
        assert_eq!(
            find_motif_indices_in_contig(&contig2, &motif3.reverse_complement()),
            vec![7, 13]
        );

        assert_eq!(find_motif_indices_in_contig(&contig2, &motif4), vec![3])
    }
}
