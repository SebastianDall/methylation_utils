pub struct Motif {
    pub sequence: String,
    pub mod_type: String,
    pub mod_position: u8,
}

impl Motif {
    pub fn new(sequence: &str, mod_type: &str, mod_position: u8) -> Self {
        Self {
            sequence: sequence.to_string(),
            mod_type: mod_type.to_string(),
            mod_position,
        }
    }

    pub fn reverse_complement(&self) -> Self {
        Self {
            sequence: (&self.sequence.chars().rev().collect::<String>()).to_string(),
            mod_type: self.mod_type.to_string(),
            mod_position: self.sequence.len() as u8 - self.mod_position - 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_creation() {
        let motif = Motif::new("GATC", "m", 1);
        assert_eq!(motif.sequence, "GATC");
        assert_eq!(motif.mod_position, 1);
    }

    #[test]
    fn test_motif_reverse_complement() {
        let motif = Motif::new("GATC", "m", 1);
        let reversed_complement_motif = motif.reverse_complement();
        assert_eq!(reversed_complement_motif.sequence, "CTAG");
        assert_eq!(reversed_complement_motif.mod_position, 2);
        assert_eq!(reversed_complement_motif.mod_type, motif.mod_type);
    }
}
