use anyhow::{bail, Context, Result};
use csv::ReaderBuilder;
use humantime::format_duration;
use indicatif::HumanDuration;
use log::info;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::Path,
    time::Instant,
};

use crate::{
    data::{pileup::PileupRecord, GenomeWorkspaceBuilder, MethylationRecord},
    data_load::load_contigs,
};

pub mod args;
pub mod processing;

pub use args::MethylationPatternArgs;
pub use processing::{
    calculate_contig_read_methylation_pattern, create_motifs, MotifMethylationDegree,
};

pub fn extract_methylation_pattern(args: MethylationPatternArgs) -> Result<()> {
    info!(
        "Running epimetheus 'methylation-pattern' with {} threads",
        &args.threads
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .context("Could not initialize threadpool")?;

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

    let mut builder = GenomeWorkspaceBuilder::new();

    let mut current_contig: Option<String> = None;
    let mut contigs_loaded = 0;
    let mut contigs_processed = 0;

    let mut pileup_records: Vec<PileupRecord> = Vec::new();
    let mut methylation_pattern_results: Vec<MotifMethylationDegree> = Vec::new();

    let mut batch_loading_duration = Instant::now();
    for row in rdr.deserialize::<PileupRecord>() {
        let pileup_rec = row?;

        if pileup_rec.n_valid_cov < args.min_valid_read_coverage {
            continue;
        }

        let contig_id = &pileup_rec.contig;

        if current_contig.as_ref() != Some(contig_id) {
            current_contig = Some(contig_id.clone());
            contigs_loaded += 1;

            if contigs_loaded > args.batches {
                let elapsed_batch_loading_duration = batch_loading_duration.elapsed();
                info!(
                    "Loading {} contigs took: {}.",
                    &args.batches,
                    format_duration(elapsed_batch_loading_duration).to_string()
                );
                let pileup_records_batch = std::mem::take(&mut pileup_records);
                let to_methylation_records_results: Result<Vec<MethylationRecord>, anyhow::Error> =
                    pool.install(|| {
                        pileup_records_batch
                            .into_par_iter()
                            .map(|pr| pr.try_into())
                            .collect()
                    });
                let mut methylation_records = to_methylation_records_results?;

                for meth_rec in methylation_records.drain(..) {
                    builder.add_record(meth_rec)?;
                }

                let workspace = builder.build();

                info!("Calculating methylation patten.");
                let calculate_methylation_pattern_duration = Instant::now();
                let mut methylation_pattern =
                    calculate_contig_read_methylation_pattern(workspace, motifs.clone(), &pool)?;
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

            let contig = match contigs.get(contig_id) {
                Some(contig) => contig,
                None => bail!("Contig not found in assembly: {contig_id}"),
            };
            builder.add_contig(contig.clone())?;
        }

        // let methylation_record = pileup_rec.try_into()?;

        // methylation_records.push(methylation_record);
        pileup_records.push(pileup_rec);
    }

    if !pileup_records.is_empty() {
        let to_methylation_records_results: Result<Vec<MethylationRecord>, anyhow::Error> = pool
            .install(|| {
                pileup_records
                    .into_par_iter()
                    .map(|pr| pr.try_into())
                    .collect()
            });
        let mut methylation_records = to_methylation_records_results?;
        for meth_rec in methylation_records.drain(..) {
            builder.add_record(meth_rec)?;
        }
        let workspace = builder.build();

        let mut methylation_pattern =
            calculate_contig_read_methylation_pattern(workspace, motifs.clone(), &pool)?;

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
