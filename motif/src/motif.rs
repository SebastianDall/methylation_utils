use crate::{IupacBase, ModType};
use anyhow::{bail, Result};

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
