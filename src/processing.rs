use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use methylome::{find_motif_indices_in_contig, motif::Motif};
use rayon::prelude::*;
use std::{
    fmt::Write,
    sync::Arc,
    str::FromStr,
};

use crate::data::{methylation::MethylationCoverage, GenomeWorkspace};

pub struct MotifMethylationDegree {
    pub contig: String,
    pub motif: Motif,
    pub median: f64,
    pub mean_read_cov: f64,
    pub n_motif_obs: u32,
    pub motif_occurences_total: u32,
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

    let results: Vec<MotifMethylationDegree> = contigs.contigs.par_iter().flat_map(|(contig_id, contig)| {
     let contig_seq = &contig.sequence;

     let mut local_results = Vec::new();

     for motif in motifs.iter() {
         let mod_type = motif.mod_type;

         let fwd_indices: Vec<usize> = find_motif_indices_in_contig(&contig_seq, motif);
         let rev_indices: Vec<usize> = find_motif_indices_in_contig(&contig_seq, &motif.reverse_complement());

         if fwd_indices.is_empty() && rev_indices.is_empty() {
             continue;
         }

         // This is the actual number of motifs in the contig
         let motif_occurences_total = fwd_indices.len() as u32 + rev_indices.len() as u32;

         let mut fwd_methylation = contig.get_methylated_positions(&fwd_indices, methylome::Strand::Positive, mod_type);
         let mut rev_methylation = contig.get_methylated_positions(&rev_indices, methylome::Strand::Negative, mod_type);

         fwd_methylation.append(&mut rev_methylation);

         let methylation_data: Vec<MethylationCoverage> = fwd_methylation.into_iter().filter_map(|maybe_cov| maybe_cov.cloned()).collect();

         if methylation_data.is_empty() {
             continue;
         }

         // This is number of motif obervations with methylation data
         let n_motif_obs = methylation_data.len() as u32;
         
         let mean_read_cov = {
             let total_cov: u64 = methylation_data.iter().map(|cov| cov.get_n_valid_cov() as u64).sum();
             total_cov as f64 / methylation_data.len() as f64
         };

         let mut fractions: Vec<f64> = methylation_data
            .iter()
            .map(|cov| cov.fraction_modified())
            .collect();

        fractions.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = if fractions.len() % 2 == 0 {
            let mid = fractions.len() / 2;
            (fractions[mid - 1] + fractions[mid]) / 2.0
        } else {
            fractions[fractions.len() / 2]
        };

         local_results.push(MotifMethylationDegree {
             contig: contig_id.clone(),
             motif: motif.clone(),
             median,
             mean_read_cov,
             n_motif_obs,
             motif_occurences_total,
         })
     }

     pb.inc(1);

     local_results

        
    }).collect();

    pb.finish_with_message("Finished processing all contigs.");

    Ok(results)
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
    use tempfile::NamedTempFile;
    use std::io::Write;

    use crate::data::contig::Contig;

    use super::*;

    #[test]
    fn test_calculate_methylation() -> Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;



        let mut workspace = GenomeWorkspace::new();

        // Add a mock contig to the workspace
        workspace.add_contig(Contig::new("contig_3".to_string(), "TGGACGATCCCGATC".to_string())).unwrap();
        let motifs = vec![
            Motif::new("GATC", "a", 1).unwrap(),
            Motif::new("GATC", "m", 3).unwrap(),
            Motif::new("GATC", "21839", 3).unwrap(),
        ];

        workspace.populate_methylation_from_pileup(pileup_file, 1).unwrap();

        let contig_methylation_pattern =
            calculate_contig_read_methylation_pattern(workspace,  motifs, 1).unwrap();

        

        let expected_median_result = vec![0.625, 1.0];
        let meth_result: Vec<f64> = contig_methylation_pattern.iter().map(|res| res.median).collect();
        assert_eq!(
            meth_result,
            expected_median_result
        );

        let expected_mean_read_cov = vec![18.75, 20.0];
        let meth_result: Vec<f64> = contig_methylation_pattern.iter().map(|res| res.mean_read_cov).collect();
        assert_eq!(
            meth_result,
            expected_mean_read_cov
        );

        let expected_n_motif_obs = vec![4, 1];
        let meth_result: Vec<u32> = contig_methylation_pattern.iter().map(|res| res.n_motif_obs).collect();
        assert_eq!(meth_result, expected_n_motif_obs);

        Ok(())
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
