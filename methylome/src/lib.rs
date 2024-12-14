use regex::Regex;

mod iupac;
pub mod modtype;
pub mod motif;

use iupac::IupacBase;
pub use modtype::ModType;
pub use motif::Motif;

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
