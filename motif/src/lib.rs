use anyhow::{bail, Result};
use regex::Regex;

mod modtype;
use modtype::ModType;

mod iupac;
use iupac::IupacBase;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Motif {
    pub sequence: Vec<IupacBase>,
    pub mod_type: ModType,
    pub mod_position: u8,
}

impl Motif {
    pub fn new(sequence: &str, mod_type: &str, mod_position: u8) -> Result<Self> {
        let mod_type = ModType::from_str(mod_type)?;

        let parsed_sequence = sequence
            .chars()
            .map(|b| {
                IupacBase::parse_char(b).map_err(|_| {
                    anyhow::anyhow!(
                        "Base '{}' in sequence '{}' is not a valid IUPAC code",
                        b,
                        sequence
                    )
                })
            })
            .collect::<Result<Vec<IupacBase>>>()?;

        if mod_position as usize > parsed_sequence.len() - 1 {
            bail!(
                "mod_position {} is out of bounds for sequence of length {}. Note mod_position is 0-indexed.",
                mod_position,
                parsed_sequence.len()
            );
        }

        let base_at_position = &parsed_sequence[mod_position as usize];
        match mod_type {
            ModType::SixMA => {
                if *base_at_position != IupacBase::A {
                    bail!(
                        "mod_position {} points to base '{}' which is invalid for 6mA.",
                        mod_position,
                        base_at_position
                    );
                }
            }
            ModType::FiveMC | ModType::FourMC => {
                if *base_at_position != IupacBase::C {
                    bail!(
                        "mod_position {} points to base '{}' which is invalid for {} modification type.",
                        mod_position, base_at_position, mod_type
                    );
                }
            }
        }

        Ok(Self {
            sequence: parsed_sequence,
            mod_type,
            mod_position,
        })
    }

    pub fn reverse_complement(&self) -> Self {
        Self {
            // sequence: (&self.sequence.chars().rev().collect::<String>()).to_string(),
            sequence: self
                .sequence
                .iter()
                .rev()
                .map(IupacBase::to_complement_base)
                .collect(),
            mod_type: self.mod_type.clone(),
            mod_position: self.sequence.len() as u8 - self.mod_position - 1,
        }
    }

    pub fn to_regex(&self) -> String {
        self.sequence.iter().map(IupacBase::to_regex).collect()
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

    fn parse_iupac_sequence(sequence: &str) -> Vec<IupacBase> {
        sequence
            .chars()
            .map(|c| IupacBase::parse_char(c).unwrap())
            .collect()
    }

    #[test]
    fn test_motif_creation() {
        let motif = Motif::new("GATC", "a", 1).unwrap();
        assert_eq!(motif.sequence, parse_iupac_sequence("GATC"));
        assert_eq!(motif.mod_type, ModType::SixMA);
        assert_eq!(motif.mod_position, 1);
    }

    #[test]
    fn test_out_of_bounds() {
        let result = Motif::new("GATC", "m", 4);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "mod_position 4 is out of bounds for sequence of length 4. Note mod_position is 0-indexed."
        );
    }

    #[test]
    fn test_unidentified_motif_type() {
        let result = Motif::new("GATC", "d", 1);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().to_string(), "Unsupported mod type: d");
    }

    #[test]
    fn test_invalid_mod_position_base() {
        let result = Motif::new("ATCG", "m", 3); // 'G' is invalid for 5mC
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "mod_position 3 points to base 'G' which is invalid for 5mC (m) modification type."
        );
    }

    #[test]
    fn test_invalid_iupac_base() {
        let result = Motif::new("ATZG", "a", 0); // 'G' is invalid for 5mC
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Base 'Z' in sequence 'ATZG' is not a valid IUPAC code"
        );
    }

    #[test]
    fn test_motif_reverse_complement() {
        let motif1 = Motif::new("GATC", "m", 3).unwrap();
        let motif2 = Motif::new("TCCCG", "m", 1).unwrap();
        let motif3 = Motif::new("RGATCY", "a", 2).unwrap();
        assert_eq!(
            motif1.reverse_complement().sequence,
            parse_iupac_sequence("GATC")
        );
        assert_eq!(
            motif2.reverse_complement().sequence,
            parse_iupac_sequence("CGGGA")
        );
        assert_eq!(
            motif3.reverse_complement().sequence,
            parse_iupac_sequence("RGATCY")
        );
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
