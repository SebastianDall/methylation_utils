use anyhow::{anyhow, bail, Context, Result};
use csv::{ReaderBuilder, StringRecord};
use humantime::format_duration;
use indicatif::HumanDuration;
use log::info;
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::Path,
    time::Instant,
};

use crate::{
    data::{GenomeWorkspaceBuilder, MethylationRecord},
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
    let mut rdr = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .flexible(false)
        .from_reader(reader);
    let mut record = StringRecord::with_capacity(100, 18);

    let mut builder = GenomeWorkspaceBuilder::new();

    let mut current_contig_id: Option<String> = None;
    let mut contigs_loaded = 0;
    let mut contigs_processed = 0;

    let mut methylation_records: Vec<MethylationRecord> = Vec::new();
    let mut methylation_pattern_results: Vec<MotifMethylationDegree> = Vec::new();

    let mut batch_loading_duration = Instant::now();
    while rdr.read_record(&mut record)? {
        let n_valid_cov: u32 = record
            .get(9)
            .ok_or_else(|| anyhow!("Missing n_valid_coverage field"))?
            .parse()
            .map_err(|_| anyhow!("Invalid coverage number."))?;
        if n_valid_cov < args.min_valid_read_coverage {
            continue;
        }

        let contig_id = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing contig field"))?
            .to_string();

        let mut current_contig_loaded = if current_contig_id.as_ref() != Some(&contig_id) {
            current_contig_id = Some(contig_id.clone());
            contigs_loaded += 1;

            contigs
                .get(&contig_id)
                .with_context(|| format!("Contig not found in assembly: {contig_id}"))?
                .clone()
            // builder.add_contig(contig.clone())?;
        } else {
            bail!("Get out")
        };

        if contigs_loaded > args.batches {
            let elapsed_batch_loading_duration = batch_loading_duration.elapsed();
            info!(
                "Loading {} contigs took: {}.",
                &args.batches,
                format_duration(elapsed_batch_loading_duration).to_string()
            );
            for meth_rec in methylation_records.drain(..) {
                builder.add_record(meth_rec)?;
            }

            let workspace = builder.build();

            info!("Calculating methylation patten.");
            let calculate_methylation_pattern_duration = Instant::now();
            let mut methylation_pattern =
                calculate_contig_read_methylation_pattern(workspace, motifs.clone(), args.threads)?;
            let elapsed_calculate_methylation_pattern_duration =
                calculate_methylation_pattern_duration.elapsed();
            info!(
                "Calculating methylation pattern took: {} - ({})",
                HumanDuration(elapsed_calculate_methylation_pattern_duration).to_string(),
                format_duration(elapsed_calculate_methylation_pattern_duration).to_string()
            );

            methylation_pattern_results.append(&mut methylation_pattern);

            contigs_processed += contigs_loaded - 1;
            info!("Finished processing {}", contigs_processed);

            builder = GenomeWorkspaceBuilder::new();
            batch_loading_duration = Instant::now();
            contigs_loaded = 1;
        }

        let methylation_record = parse_to_methylation_record(contig_id, n_valid_cov, &record)?;
        current_contig_loaded.add_methylation_record(methylation_record)?;
    }

    if !methylation_records.is_empty() {
        for meth_rec in methylation_records.drain(..) {
            builder.add_record(meth_rec)?;
        }
        let workspace = builder.build();

        let mut methylation_pattern =
            calculate_contig_read_methylation_pattern(workspace, motifs.clone(), args.threads)?;

        methylation_pattern_results.append(&mut methylation_pattern);
        contigs_processed += contigs_loaded;
        info!("Finished loading {} contigs", contigs_processed);
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
