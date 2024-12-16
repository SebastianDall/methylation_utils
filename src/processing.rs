use anyhow::{Context, Result};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressState, ProgressStyle};
use log::{error, info};
use methylome::{find_motif_indices_in_contig, motif::Motif};
use polars::{datatypes::DataType, frame::DataFrame, lazy::frame::LazyFrame, prelude::*};
use rayon::prelude::*;
use std::{
    env,
    fmt::Write,
    sync::{Arc, Mutex},
    time::Duration,
    str::FromStr,
};

use crate::data::GenomeWorkspace;

struct MotifMethylationDegree {
    contig: String,
    motif: Motif,
    median: f64,
    mean_read_coverage: f64,
}

pub fn calculate_contig_read_methylation_pattern(
    contigs: GenomeWorkspace,
    motifs: Vec<Motif>,
    num_threads: usize,
) -> Result<Vec<MotifMethylationDegree>> {

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Could not initialize threadpool");

    let tasks = contigs.contigs.len() as u64;
    let pb = ProgressBar::new(tasks);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
        )
        .unwrap()
        .with_key("eta", |state: &ProgressState, w: &mut dyn Write| {
            write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap()
        })
        .progress_chars("#>-"),
    );

    let motifs = Arc::new(motifs);

    let results: <Vec<MotifMethylationDegree>> = contigs.contigs.par_iter().flat_map()
     let contig_seq = &contig.sequence;

     let mut local_results = Vec::new();

     for motif in motifs.iter() {
         let mod_type = motif.mod_type;

         let fwd_indices = find_motif_indices_in_contig(&contig_seq, motif);
         let rev_indices = find_motif_indices_in_contig(&contig_seq, &motif.reverse_complement());

         if fwd_indices.is_empty() && rev_indices.is_empty() {
             continue;
         }

         

         
         
     }
        
    });

    

    Ok()
}

pub fn create_motifs(motifs_str: Vec<String>) -> Result<Vec<Motif>> {
    motifs_str.into_iter().map(|motif| {
        let parts: Vec<&str> = motif.split("_").collect();

        if parts.len() != 3 {
            anyhow::bail!(
                "Invalid motif format '{}' encountered. Expected format: '<sequence>_<mod_type>_<mod_position>'",
                motif
            );
        }

            let sequence = parts[0];
            let mod_type = parts[1];
            let mod_position = u8::from_str(parts[2]).with_context(|| {
                format!("Failed to parse mod_position '{}' in motif '{}'.", parts[2], motif)
            })?;

            Motif::new(sequence, mod_type, mod_position).with_context(|| {
                format!("Failed to create motif from '{}'", motif)
            })
        
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_methylation() {
        let subpileup = df!(
            "contig" => ["contig_3","contig_3","contig_3","contig_3","contig_3"],
            "strand" => ["+", "+", "+", "-", "-"],
            "mod_type" => ["a", "m", "a", "a", "a"],
            "start" => [6, 8, 12, 7, 13],
            "N_modified" => [15, 20, 5, 20 ,5],
            "N_valid_cov" => [15, 20, 20, 20, 20],
            "motif_mean" => [1.0, 1.0, 0.25, 1.0, 0.25]
        )
        .expect("Could not initialize dataframe");

        let subpileups = vec![subpileup];

        let mut contig_map = ContigMap::new();
        contig_map.insert("contig_3".to_string(), "TGGACGATCCCGATC".to_string());

        let motifs = vec![
            Motif::new("GATC", "a", 1).unwrap(),
            Motif::new("GATC", "m", 3).unwrap(),
            Motif::new("GATC", "21839", 3).unwrap(),
        ];

        let contig_methylation_pattern =
            calculate_contig_read_methylation_pattern(contig_map, subpileups, motifs, 1).unwrap();

        println!("{:#?}", contig_methylation_pattern);

        let expected_result = Column::new("median".into(), [0.625, 1.0]);
        assert_eq!(
            contig_methylation_pattern.column("median").unwrap(),
            &expected_result
        );

        let expected_mean_read_cov = Column::new("mean_read_cov".into(), [18.75, 20.0]);
        assert_eq!(
            contig_methylation_pattern.column("mean_read_cov").unwrap(),
            &expected_mean_read_cov
        )
    }

    #[test]
    fn test_create_motifs_success() {
        let motifs_args = vec!["GATC_a_1".to_string()];
        let result = create_motifs(motifs_args);
        assert!(result.is_ok(), "Expected Ok, but got err: {:?}", result.err());
    }
    #[test]
    fn test_create_motifs_failure() {
        let motifs_args = vec!["GATC_a_3".to_string()];
        let result = create_motifs(motifs_args);
        assert!(result.is_err(), "Expected Err, but got Ok: {:?}", result.ok());
    }
    
}
