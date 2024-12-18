use anyhow::{Context, Result};
use clap::Parser;
use humantime::format_duration;
use indicatif::HumanDuration;
use log::info;
use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
    time::Instant,
};

mod argparser;
mod data;
mod data_load;
mod processing;

use argparser::Args;
use data_load::load_contigs;
use processing::{calculate_contig_read_methylation_pattern, create_motifs};

fn main() -> Result<()> {
    // let guard = pprof::ProfilerGuard::new(1000).unwrap();
    let total_duration = Instant::now();
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();
    info!("Running methylation_utils with {} threads", &args.threads);

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

    info!("Loading assembly");
    let mut contigs = load_contigs(&args.assembly)
        .with_context(|| format!("Error loading assembly from path: '{}'", args.assembly))?;

    if contigs.contigs.len() == 0 {
        anyhow::bail!("No contigs are loaded!");
    }

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

    info!("Loading methylation from pileup.");
    contigs.populate_methylation_from_pileup(args.pileup, args.min_valid_read_coverage)?;

    info!("Removing contigs with no methylation");
    contigs.prune_empty_contigs();

    let elapsed_preparation_time = preparation_duration.elapsed();
    info!(
        "Data preparation took: {} - ({})",
        HumanDuration(elapsed_preparation_time).to_string(),
        format_duration(elapsed_preparation_time).to_string()
    );

    info!("Finding contig methylation pattern");
    let finding_methylation_pattern_duration = Instant::now();
    let mut contig_methylation_pattern =
        calculate_contig_read_methylation_pattern(contigs, motifs, args.threads)
            .context("Unable to find methyation pattern.")?;

    let elapsed_finding_methylation_pattern_duration =
        finding_methylation_pattern_duration.elapsed();
    info!(
        "Methylation pattern took: {} - ({})",
        HumanDuration(elapsed_finding_methylation_pattern_duration).to_string(),
        format_duration(elapsed_finding_methylation_pattern_duration).to_string()
    );

    contig_methylation_pattern.sort_by(|a, b| a.contig.cmp(&b.contig));

    let outfile = std::fs::File::create(outpath)
        .with_context(|| format!("Failed to create file at: {:?}", outpath))?;
    let mut writer = BufWriter::new(outfile);

    writeln!(
        writer,
        "contig\tmotif\tmod_type\tmod_position\tmedian\tmean_read_cov\tN_motif_obs\tmotif_occurences_total"
    )?;

    for entry in &contig_methylation_pattern {
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
