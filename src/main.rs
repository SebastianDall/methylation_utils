use clap::Parser;
use core::panic;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use motif::{find_motif_indices_in_contig, Motif};
use polars::{
    datatypes::DataType,
    error::PolarsResult,
    frame::DataFrame,
    lazy::frame::{LazyCsvReader, LazyFileListReader, LazyFrame},
    prelude::*,
};
use rayon::prelude::*;
use seq_io::fasta::{Reader, Record};
use std::{
    collections::HashMap,
    fmt::Write,
    fs,
    path::Path,
    process, str,
    sync::{Arc, Mutex},
};

type ContigMap = HashMap<String, String>;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long, required = true)]
    pileup: String,

    #[arg(short, long, required = true)]
    assembly: String,

    #[arg(short, long, required = true)]
    output: String,

    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    #[arg(short, long, required = true, num_args(1..))]
    motifs: Option<Vec<String>>,

    #[arg(long, default_value_t = 3)]
    min_valid_read_coverage: u8,
}

fn load_pileup_lazy<P: AsRef<Path>>(path: P) -> PolarsResult<LazyFrame> {
    let old_column_names: Vec<String> = (1..19).map(|c| format!("column_{}", c)).collect();
    let new_column_names = vec![
        "contig",
        "start",
        "end",
        "mod_type",
        "score",
        "strand",
        "start2",
        "end2",
        "color",
        "N_valid_cov",
        "percent_modified",
        "N_modified",
        "N_canonical",
        "N_other_mod",
        "N_delete",
        "N_fail",
        "N_diff",
        "N_nocall",
    ];

    let lf_pileup = LazyCsvReader::new(path)
        .with_has_header(false)
        .with_separator(b'\t')
        .finish()?
        .rename(&old_column_names, &new_column_names, true);

    Ok(lf_pileup)
}

fn load_contigs<P: AsRef<Path>>(path: P) -> Result<ContigMap, Box<dyn std::error::Error>> {
    let mut fasta_reader = Reader::from_path(path).unwrap();
    let mut contigs = ContigMap::new();

    while let Some(record_result) = fasta_reader.next() {
        let record = record_result?;

        let id = record.id()?.to_string();
        let seq = str::from_utf8(record.seq())?.to_string();

        contigs.insert(id, seq);
    }
    Ok(contigs)
}

fn create_motifs(motifs: Vec<String>) -> Vec<Motif> {
    let mut motifs_as_struct = Vec::new();

    for motif in motifs {
        let parts: Vec<&str> = motif.split("_").collect();

        if parts.len() == 3 {
            let sequence = parts[0];
            let mod_type = parts[1];
            let mod_position: u8 = parts[2].parse().unwrap();

            let motif = Motif::new(sequence, mod_type, mod_position);
            motifs_as_struct.push(motif);
        }
    }
    motifs_as_struct
}

fn calculate_contig_read_methylation_pattern(
    contigs: ContigMap,
    pileup: DataFrame,
    motifs: Vec<Motif>,
    num_threads: usize,
) -> DataFrame {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Failed to build thread pool");

    let tasks = contigs.keys().cloned().collect::<Vec<String>>().len() as u64;
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

    let pb_arcmut = Arc::new(Mutex::new(pb));

    let motifs = Arc::new(motifs);
    let pileup = Arc::new(pileup.lazy());

    let results: Vec<LazyFrame> = contigs
        .par_iter()
        .map(|(contig_id, contig_seq)| {
            let subpileup = pileup
                .as_ref()
                .clone()
                .filter(col("contig").eq(lit(contig_id.clone())))
                .collect()
                .unwrap();

            // let mut local_read_methylation_df = DataFrame::empty();
            let f1: Field = Field::new("contig".into(), DataType::String);
            let f2: Field = Field::new("median".into(), DataType::Float64);
            let f3: Field = Field::new("N_motif_obs".into(), DataType::UInt32);
            let f4: Field = Field::new("motif".into(), DataType::String);
            let f5: Field = Field::new("mod_type".into(), DataType::String);
            let f6: Field = Field::new("mod_position".into(), DataType::Int32);

            let schema = Schema::from_iter(vec![f1, f2, f3, f4, f5, f6]);
            let mut local_read_methylation_df = DataFrame::empty_with_schema(&schema);

            for motif in motifs.iter() {
                let fwd_indices = find_motif_indices_in_contig(&contig_seq, &motif);

                let rev_indices =
                    find_motif_indices_in_contig(&contig_seq, &motif.reverse_complement());

                let p_fwd = subpileup
                    .clone()
                    .lazy()
                    .filter(
                        col("strand")
                            .eq(lit("+"))
                            .and(col("mod_type").eq(lit(motif.mod_type.clone())))
                            .and(
                                col("start").is_in(lit(Series::new("start".into(), &fwd_indices))),
                            ),
                    )
                    .collect()
                    .unwrap();

                let p_rev = subpileup
                    .clone()
                    .lazy()
                    .filter(
                        col("strand")
                            .eq(lit("-"))
                            .and(col("mod_type").eq(lit(motif.mod_type.clone())))
                            .and(
                                col("start").is_in(lit(Series::new("start".into(), &rev_indices))),
                            ),
                    )
                    .collect()
                    .unwrap();

                let p_con = concat(
                    [p_fwd.clone().lazy(), p_rev.clone().lazy()],
                    UnionArgs::default(),
                )
                .unwrap()
                .collect()
                .unwrap();

                let p_read_methylation = p_con
                    .lazy()
                    .with_columns([(col("N_modified") / col("N_valid_cov")).alias("motif_mean")])
                    .group_by([col("contig")])
                    .agg([
                        (col("motif_mean").median()).alias("median"),
                        (col("contig").count()).alias("N_motif_obs"),
                    ])
                    .with_columns([
                        (lit(motif.sequence.clone())).alias("motif"),
                        (lit(motif.mod_type.clone())).alias("mod_type"),
                        (lit(motif.mod_position.clone() as i8)).alias("mod_position"),
                    ])
                    .collect()
                    .unwrap();
                local_read_methylation_df = concat(
                    [local_read_methylation_df.lazy(), p_read_methylation.lazy()],
                    UnionArgs::default(),
                )
                .unwrap()
                .collect()
                .unwrap();
            }

            {
                let pb = pb_arcmut.lock().unwrap();
                pb.inc(1);
            }

            local_read_methylation_df.lazy()
        })
        .collect();

    let read_methylation_df = concat(&results, UnionArgs::default())
        .unwrap()
        .collect()
        .unwrap();

    pb_arcmut
        .lock()
        .unwrap()
        .finish_with_message("Finished processing contigs");

    read_methylation_df
}

fn main() {
    // let guard = pprof::ProfilerGuard::new(100).unwrap();
    let args = Args::parse();

    let outpath = Path::new(&args.output);

    match outpath.extension() {
        Some(ext) if ext == "tsv" => {
            if let Some(parent) = outpath.parent() {
                fs::create_dir_all(parent).expect("Cannot create output dir");
            }
        }
        Some(ext) => {
            println!("Incorrect file extension: {:#?}. Should be tsv", ext);
            process::exit(1);
        }
        None => {
            println!("No filename provided");
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

    let lf_pileup = load_pileup_lazy(&args.pileup).expect("Error loading pileup");

    let contigs = load_contigs(&args.assembly).expect("Error loading assembly");
    let contig_ids: Vec<String> = contigs.keys().cloned().collect();

    let pileup = lf_pileup
        .select([
            col("contig"),
            col("mod_type"),
            col("strand"),
            col("start"),
            col("N_valid_cov"),
            col("N_modified"),
        ])
        .filter(
            col("N_valid_cov")
                .gt_eq(lit(args.min_valid_read_coverage as i64))
                .and(col("contig").is_in(lit(Series::new("contig".into(), &contig_ids)))),
        )
        .collect()
        .expect("Error collecting pileup");

    let mut contig_methylation_pattern =
        calculate_contig_read_methylation_pattern(contigs, pileup, motifs, args.threads);

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
    //     use std::fs::File;
    //     use std::io::Write;

    //     let mut file = File::create("flamegraph.svg").unwrap();
    //     report.flamegraph(&mut file).unwrap();
    // }
}
