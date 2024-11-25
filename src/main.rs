use clap::Parser;
use core::panic;
use polars::prelude::*;
use std::{env, fs, path::Path, process};

mod data_load;
use data_load::{load_contigs, load_pileup_lazy};

mod types;

mod argparser;
use argparser::Args;

mod processing;
use processing::{calculate_contig_read_methylation_pattern, create_motifs, create_subpileups};

fn main() {
    // let guard = pprof::ProfilerGuard::new(100).unwrap();
    let args = Args::parse();
    println!("Running methylation_utils with {} threads", &args.threads);

    env::set_var("POLARS_MAX_THREADS", &args.threads.to_string());

    match env::var("POLARS_MAX_THREADS") {
        Ok(_val) => {}
        Err(e) => println!("POLARS_MAX_THREADS is not set: {}", e),
    }

    let outpath = Path::new(&args.output);

    match outpath.extension() {
        Some(ext) if ext == "tsv" => {
            if let Some(parent) = outpath.parent() {
                fs::create_dir_all(parent).expect("Cannot create output dir");
            }
        }
        Some(ext) => {
            eprintln!("Incorrect file extension: {:#?}. Should be tsv", ext);
            process::exit(1);
        }
        None => {
            eprintln!("No filename provided");
            process::exit(1);
        }
    }

    let motifs = match args.motifs {
        Some(motifs) => {
            println!("Motifs loaded");
            motifs
        }
        _ => panic!("No motifs found"),
    };

    let motifs = create_motifs(motifs);

    println!("Loading pileup");
    let lf_pileup = load_pileup_lazy(&args.pileup).expect("Error loading pileup");

    println!("Loading assembly");
    let contigs = load_contigs(&args.assembly).expect("Error loading assembly");
    let contig_ids: Vec<String> = contigs.keys().cloned().collect();

    let subpileups = create_subpileups(lf_pileup, contig_ids, args.min_valid_read_coverage);

    println!("Finding contig methylation pattern");
    let mut contig_methylation_pattern =
        calculate_contig_read_methylation_pattern(contigs, subpileups, motifs, args.threads);

    match std::fs::File::create(outpath) {
        Ok(mut f) => {
            CsvWriter::new(&mut f)
                .include_header(true)
                .with_separator(b'\t')
                .finish(&mut contig_methylation_pattern)
                .unwrap();
        }
        Err(e) => {
            println!("Error writing tsv file: {:?}", e);
        }
    };
    // if let Ok(report) = guard.report().build() {
    // use std::fs::File;
    //     use std::io::Write;

    //     let mut file = File::create("flamegraph.svg").unwrap();
    //     report.flamegraph(&mut file).unwrap();
    // }
}
