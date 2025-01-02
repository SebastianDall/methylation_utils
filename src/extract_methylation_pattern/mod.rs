use anyhow::{Context, Result};
use batch_loader::BatchLoader;
use humantime::format_duration;
use indicatif::HumanDuration;
use log::{debug, error, info};
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::Path,
    time::Instant,
};

use crate::{
    data_load::load_contigs,
    processing::{
        calculate_contig_read_methylation_pattern, create_motifs, MotifMethylationDegree,
    },
};

pub mod args;
pub mod batch_loader;
pub mod utils;

pub use args::MethylationPatternArgs;
pub use utils::parse_to_methylation_record;

pub fn extract_methylation_pattern(args: MethylationPatternArgs) -> Result<()> {
    info!(
        "Running epimetheus 'methylation-pattern' with {} threads",
        &args.threads
    );

    let outpath = Path::new(&args.output);

    if let Some(ext) = outpath.extension() {
        if ext != "tsv" {
            anyhow::bail!("Incorrect file extension {:?}. Should be tsv", ext);
        }
        if let Some(parent) = outpath.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Could not create parent directory: {:?}", parent))?;
        }
    } else {
        anyhow::bail!("No filename provided for output. Should be a .tsv file.");
    }

    let motifs = match args.motifs {
        Some(motifs) => {
            info!("Motifs loaded");
            motifs
        }
        _ => {
            anyhow::bail!("No motifs found");
        }
    };

    let motifs = create_motifs(motifs).context("Failed to parse motifs")?;
    info!("Successfully parsed motifs.");

    info!("Loading assembly");
    let contigs = load_contigs(&args.assembly)
        .with_context(|| format!("Error loading assembly from path: '{}'", args.assembly))?;

    if contigs.len() == 0 {
        anyhow::bail!("No contigs are loaded!");
    }
    info!("Total contigs in assembly: {}", contigs.len());

    info!("Processing Pileup");
    let file = File::open(&args.pileup)?;
    let reader = BufReader::new(file);

    let batch_loader = BatchLoader::new(
        reader,
        contigs,
        args.batch_size,
        args.min_valid_read_coverage,
    );

    let mut methylation_pattern_results: Vec<MotifMethylationDegree> = Vec::new();

    let mut batch_processing_time = Instant::now();
    let mut contigs_processed = 0;
    for ws_result in batch_loader {
        match ws_result {
            Ok(workspace) => {
                debug!("Workspace initialized");
                let contigs_in_batch = workspace.get_workspace().len() as u32;
                let mut methylation_pattern = calculate_contig_read_methylation_pattern(
                    workspace,
                    motifs.clone(),
                    args.threads,
                )?;
                methylation_pattern_results.append(&mut methylation_pattern);

                contigs_processed += contigs_in_batch;
                let elapsed_batch_processing_time = batch_processing_time.elapsed();
                info!(
                    "Finished processing {} contigs. Processing time: {}",
                    contigs_processed,
                    format_duration(elapsed_batch_processing_time).to_string()
                );
                batch_processing_time = Instant::now();
            }
            Err(e) => {
                error!("Error reading batch: {e}");
            }
        }
    }

    methylation_pattern_results.sort_by(|a, b| a.contig.cmp(&b.contig));

    let outfile = std::fs::File::create(outpath)
        .with_context(|| format!("Failed to create file at: {:?}", outpath))?;
    let mut writer = BufWriter::new(outfile);

    writeln!(
        writer,
        "contig\tmotif\tmod_type\tmod_position\tmedian\tmean_read_cov\tN_motif_obs\tmotif_occurences_total"
    )?;

    for entry in &methylation_pattern_results {
        let motif_sequence = entry.motif.sequence_to_string();
        let mod_type_str = entry.motif.mod_type.to_pileup_code();
        let mod_position = entry.motif.mod_position;

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            entry.contig,
            motif_sequence,
            mod_type_str,
            mod_position,
            entry.median,
            entry.mean_read_cov,
            entry.n_motif_obs,
            entry.motif_occurences_total
        )?;

        writer.flush()?;
    }

    Ok(())
}
