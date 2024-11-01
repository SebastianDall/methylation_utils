use phf::phf_map;
use regex::Regex;

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
            mod_type: self.mod_type.to_string(),
            mod_position: self.sequence.len() as u8 - self.mod_position - 1,
        }
    }

    pub fn to_regex(&self) -> String {
        self.sequence
            .chars()
            .map(|c| match IUPAC_MAPPING.get(&c.to_string().as_str()) {
                Some(s) => {
                    format!("[{}]", s)
                }
                None => c.to_string(),
            })
            .collect()
    }
}

pub fn find_motif_indices_in_contig(contig: &str, motif: &Motif) -> Vec<u32> {
    let regex_str = motif.to_regex();
    let re = Regex::new(&regex_str).expect("Expected regex patter");

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
        let motif = Motif::new("GATC", "m", 1);
        assert_eq!(motif.sequence, "GATC");
        assert_eq!(motif.mod_position, 1);
    }

    #[test]
    fn test_motif_reverse_complement() {
        let motif1 = Motif::new("GATC", "m", 3);
        let motif2 = Motif::new("TCCCG", "m", 1);
        let motif3 = Motif::new("RGATCY", "a", 2);
        assert_eq!(motif1.reverse_complement().sequence, "GATC");
        assert_eq!(motif2.reverse_complement().sequence, "CGGGA");
        assert_eq!(motif3.reverse_complement().sequence, "RGATCY");
        assert_eq!(motif1.reverse_complement().mod_type, "m");
        assert_eq!(motif2.reverse_complement().mod_type, "m");
        assert_eq!(motif3.reverse_complement().mod_type, "a");
        assert_eq!(motif1.reverse_complement().mod_position, 0);
        assert_eq!(motif2.reverse_complement().mod_position, 3);
        assert_eq!(motif3.reverse_complement().mod_position, 3);
    }

    #[test]
    fn test_to_regex() {
        let motif1 = Motif::new("GATC", "m", 3);
        let motif2 = Motif::new("RGATCY", "m", 4);

        assert_eq!(motif1.to_regex(), "GATC");
        assert_eq!(motif2.to_regex(), "[AG]GATC[CT]");
    }

    #[test]
    fn test_find_motif_indices_in_contig() {
        let contig = "GGATCTCCATGATC".to_string();
        let contig2 = "TGGACGATCCCGATC".to_string();
        let motif1 = Motif::new("GATC", "m", 3);
        let motif2 = Motif::new("RGATCY", "m", 4);
        let motif3 = Motif::new("GATC", "a", 1);

        assert_eq!(find_motif_indices_in_contig(&contig, &motif1), vec![4, 13]);
        assert_eq!(find_motif_indices_in_contig(&contig, &motif2), vec![4]);

        assert_eq!(find_motif_indices_in_contig(&contig2, &motif3), vec![6, 12]);
        assert_eq!(
            find_motif_indices_in_contig(&contig2, &motif3.reverse_complement()),
            vec![7, 13]
        );
    }
}
