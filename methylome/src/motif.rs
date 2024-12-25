use crate::{IupacBase, ModType};
use anyhow::{bail, Result};
use std::str::FromStr;

/// Represents a biological motif, which includes a nucleotide sequence,
/// its modification type, and the position of the modification.
///
/// # Fields
/// - `sequence`: A vector of IUPAC bases representing the motif sequence.
/// - `mod_type`: The type of modification (e.g., 6mA, 5mC).
/// - `mod_position`: The position of the modification within the sequence (0-indexed).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Motif {
    pub sequence: Vec<IupacBase>,
    pub mod_type: ModType,
    pub mod_position: u8,
}

impl Motif {
    /// Constructs a new `Motif` from a string sequence, modification type, and modification position.
    ///
    /// # Arguments
    /// - `sequence`: A string representing the nucleotide sequence (using IUPAC codes).
    /// - `mod_type`: A string representing the modification type (e.g., "a" (6mA), "m" (5mC), "21839" (4mC)0).
    /// - `mod_position`: The 0-indexed position of the modification in the sequence.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The `sequence` contains invalid IUPAC codes.
    /// - The `mod_position` is out of bounds for the sequence.
    /// - The `mod_type` does not match the base at `mod_position` (e.g., 6mA must modify an 'A').
    ///
    /// # Examples
    /// ```
    /// use methylome::{Motif, ModType};
    ///
    /// let motif = Motif::new("GATC", "a", 1).unwrap();
    /// assert_eq!(motif.mod_type, ModType::SixMA);
    /// ```
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

    /// Returns the reverse complement of the motif.
    ///
    /// The reverse complement reverses the sequence and replaces each base
    /// with its complement (e.g., A ↔ T, C ↔ G). The modification position
    /// is adjusted to reflect its position in the reverse-complemented sequence.
    ///
    /// # Examples
    /// ```
    /// use methylome::Motif;
    ///
    /// let motif = Motif::new("TCCCG", "m", 1).unwrap();
    /// let rev_comp = motif.reverse_complement();
    /// assert_eq!(rev_comp.sequence_to_string(), "CGGGA");
    /// assert_eq!(rev_comp.mod_position, 3);
    /// ```
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

    /// Converts the motif sequence into a regular expression string.
    ///
    /// Each base in the sequence is mapped to its corresponding regex
    /// pattern based on IUPAC codes. For example, `R` (purine) becomes `[AG]`.
    ///
    /// # Examples
    /// ```
    /// use methylome::Motif;
    ///
    /// let motif = Motif::new("RGATCY", "a", 2).unwrap();
    /// let regex = motif.to_regex();
    /// assert_eq!(regex, "[AG]GATC[CT]");
    /// ```
    pub fn to_regex(&self) -> String {
        self.sequence.iter().map(IupacBase::to_regex).collect()
    }

    /// Converts the motif sequence into a plain string representation.
    ///
    /// This method maps each IUPAC base in the sequence to its corresponding character.
    ///
    /// # Examples
    /// ```
    /// use methylome::Motif;
    ///
    /// let motif = Motif::new("GATC", "m", 3).unwrap();
    /// let sequence = motif.sequence_to_string();
    /// assert_eq!(sequence, "GATC");
    /// ```
    pub fn sequence_to_string(&self) -> String {
        self.sequence.iter().map(IupacBase::to_string).collect()
    }
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
}
