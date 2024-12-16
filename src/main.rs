use anyhow::{Context, Result};
use clap::Parser;
use humantime::format_duration;
use indicatif::HumanDuration;
use log::{error, info};
use polars::prelude::*;
use std::{env, fs, path::Path, time::Instant};

mod argparser;
mod data;
mod data_load;
mod processing;
mod types;

use argparser::Args;
use data_load::{load_contigs, load_pileup_lazy};
use processing::{calculate_contig_read_methylation_pattern, create_motifs};

fn main() -> Result<()> {
    // let guard = pprof::ProfilerGuard::new(1000).unwrap();
    let total_duration = Instant::now();
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();
    info!("Running methylation_utils with {} threads", &args.threads);

    env::set_var("POLARS_MAX_THREADS", &args.threads.to_string());

    if env::var("POLARS_MAX_THREADS").is_err() {
        error!("WARNING: POLARS_MAX_THREADS is not set. This will hurt performance");
    }

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

    let preparation_duration = Instant::now();
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

    info!("Loading pileup");
    let lf_pileup = load_pileup_lazy(&args.pileup).context("Error loading pileup")?;

    info!("Loading assembly");
    let contigs = load_contigs(&args.assembly)
        .with_context(|| format!("Error loading assembly from path: '{}'", args.assembly))?;
    // let contig_ids: Vec<String> = contigs.keys().cloned().collect();

    // if contig_ids.len() == 0 {
    //     anyhow::bail!("No contigs are loaded!");
    // }

    // let batches = if args.batches == 0 {
    //     info!("Loading pileup without batching");
    //     contig_ids.len()
    // } else {
    //     info!("Loading pileup with {} contigs at a time", args.batches);
    //     args.batches
    // };

    // let subpileups =
    //     create_subpileups(lf_pileup, contig_ids, args.min_valid_read_coverage, batches)
    //         .context("Unable to create subpileups.")?;

    let elapsed_preparation_time = preparation_duration.elapsed();
    info!(
        "Data preparation took: {} - ({})",
        HumanDuration(elapsed_preparation_time).to_string(),
        format_duration(elapsed_preparation_time).to_string()
    );

    info!("Finding contig methylation pattern");
    let finding_methylation_pattern_duration = Instant::now();
    // let mut contig_methylation_pattern =
    // calculate_contig_read_methylation_pattern(contigs, subpileups, motifs, args.threads)
    //     .context("Unable to find methyation pattern.")?;
    let elapsed_finding_methylation_pattern_duration =
        finding_methylation_pattern_duration.elapsed();
    info!(
        "Methylation pattern took: {} - ({})",
        HumanDuration(elapsed_finding_methylation_pattern_duration).to_string(),
        format_duration(elapsed_finding_methylation_pattern_duration).to_string()
    );

    // let mut outfile = std::fs::File::create(outpath)
    //     .with_context(|| format!("Failed to create file at: {:?}", outpath))?;

    // CsvWriter::new(&mut outfile)
    //     .include_header(true)
    //     .with_separator(b'\t')
    //     .finish(&mut contig_methylation_pattern)
    //     .with_context(|| format!("Error writing tsv"))?;

    let elapsed_total_duration = total_duration.elapsed();
    info!(
        "Total time: {} - ({})",
        HumanDuration(elapsed_total_duration).to_string(),
        format_duration(elapsed_total_duration).to_string()
    );
    // if let Ok(report) = guard.report().build() {
    //     use std::fs::File;
    //     use std::io::Write;

    //     let mut file = File::create("flamegraph.svg").unwrap();
    //     report.flamegraph(&mut file).unwrap();
    // }
    Ok(())
}
