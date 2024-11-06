use clap::Parser;
use core::panic;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use motif::{find_motif_indices_in_contig, Motif};
use polars::{datatypes::DataType, frame::DataFrame, lazy::frame::LazyFrame, prelude::*};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fmt::Write,
    fs,
    path::Path,
    process,
    sync::{Arc, Mutex},
};

mod data_load;
use data_load::{create_motifs, load_contigs, load_pileup_lazy};

mod types;
use types::ContigMap;

mod argparser;
use argparser::Args;

fn create_subpileups(
    pileup: LazyFrame,
    contig_ids: Vec<String>,
    min_valid_read_coverage: u32,
) -> HashMap<String, DataFrame> {
    println!("Filtering pileup");
    let pileup = pileup
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
                .gt_eq(lit(min_valid_read_coverage as i64))
                .and(col("contig").is_in(lit(Series::new("contig".into(), &contig_ids)))),
        )
        .collect()
        .expect("Error collecting pileup");

    let subpileups = pileup
        .partition_by(["contig"], true)
        .expect("Couldn't partition pileup");

    drop(pileup);

    let mut subpileups_map = HashMap::new();
    for df in subpileups {
        if let Ok(contig_series) = df.column("contig") {
            if let Ok(contig_id_anyvalue) = contig_series.get(0) {
                match contig_id_anyvalue {
                    AnyValue::String(contig_id_str) => {
                        subpileups_map.insert(contig_id_str.to_string(), df);
                    }
                    AnyValue::StringOwned(contig_id_str) => {
                        subpileups_map.insert(contig_id_str.to_string(), df);
                    }
                    _ => {
                        eprintln!(
                            "Expected String value for contig ID, found {:?}",
                            contig_id_anyvalue
                        );
                    }
                }
            } else {
                eprintln!("Failed to get value from contig series");
            }
        } else {
            eprintln!("Failed to get contig column");
        }
    }
    println!(
        "Number of subpileups generated {:?}",
        subpileups_map.keys().count()
    );
    subpileups_map
}

fn calculate_contig_read_methylation_pattern(
    contigs: ContigMap,
    subpileups_map: HashMap<String, DataFrame>,
    motifs: Vec<Motif>,
    num_threads: usize,
) -> DataFrame {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Failed to build thread pool");

    let contig_ids = subpileups_map.keys().cloned().collect::<Vec<String>>();

    let tasks = contig_ids.len() as u64;
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

    let f1: Field = Field::new("contig".into(), DataType::String);
    let f2: Field = Field::new("median".into(), DataType::Float64);
    let f3: Field = Field::new("N_motif_obs".into(), DataType::UInt32);
    let f4: Field = Field::new("motif".into(), DataType::String);
    let f5: Field = Field::new("mod_type".into(), DataType::String);
    let f6: Field = Field::new("mod_position".into(), DataType::Int32);

    let schema = Schema::from_iter(vec![f1, f2, f3, f4, f5, f6]);
    let empty_df = DataFrame::empty_with_schema(&schema);

    let results = Arc::new(Mutex::new(Vec::new()));

    contig_ids.chunks(1000).for_each(|batch| {
        batch.par_iter().for_each(|contig_id| {
            let contig_seq = contigs
                .get(contig_id)
                .expect("Could not find contig in ContigMap");
            let subpileup = subpileups_map
                .get(contig_id)
                .expect("Could not find subpileup in subpileups_map");

            let mut local_read_methylation_df = empty_df.clone();

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

                if (p_fwd.shape().0, p_rev.shape().0) == (0, 0) {
                    // TODO: put in a log file.
                    // println!(
                    //     "{} not found in contig {}",
                    //     motif.sequence.clone(),
                    //     contig_id
                    // );
                    continue;
                }

                let p_con =
                    p_fwd
                        .vstack(&p_rev)
                        .unwrap()
                        .lazy()
                        .with_columns([(col("N_modified").cast(DataType::Float64)
                            / col("N_valid_cov").cast(DataType::Float64))
                        .alias("motif_mean")]);

                let p_read_methylation = p_con
                    .group_by([col("contig")])
                    .agg([
                        (col("motif_mean").median()).alias("median"),
                        (col("contig").count()).alias("N_motif_obs"),
                    ])
                    .with_columns([
                        (lit(motif.sequence.clone())).alias("motif"),
                        (lit(motif.mod_type.clone())).alias("mod_type"),
                        (lit(motif.mod_position.clone() as i8)).alias("mod_position"),
                    ]);

                local_read_methylation_df = concat(
                    [local_read_methylation_df.lazy(), p_read_methylation],
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

            let mut results_lock = results.lock().expect("Could not open results mutex");
            results_lock.push(local_read_methylation_df.lazy());
        });
    });

    let final_results: Vec<LazyFrame> = results.lock().unwrap().clone();

    let read_methylation_df = concat(&final_results, UnionArgs::default())
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
    //     use std::fs::File;
    //     use std::io::Write;

    //     let mut file = File::create("flamegraph.svg").unwrap();
    //     report.flamegraph(&mut file).unwrap();
    // }
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
            "N_modified" => [20, 20, 5, 20 ,5],
            "N_valid_cov" => [20 , 20, 20, 20, 20]
        )
        .expect("Could not initialize dataframe");

        let mut subpileups = HashMap::new();
        subpileups.insert("contig_3".to_string(), subpileup);

        let mut contig_map = ContigMap::new();
        contig_map.insert("contig_3".to_string(), "TGGACGATCCCGATC".to_string());

        let motifs = vec![Motif::new("GATC", "a", 1)];

        let contig_methylation_pattern =
            calculate_contig_read_methylation_pattern(contig_map, subpileups, motifs, 1);

        println!("{:#?}", contig_methylation_pattern);

        let expected_result = Column::new("median".into(), [0.625]);
        assert_eq!(
            contig_methylation_pattern.column("median").unwrap(),
            &expected_result
        );
    }
}
