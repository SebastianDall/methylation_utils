use anyhow::{anyhow, Context, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};
use log::{error, info};
use motif::{find_motif_indices_in_contig, Motif};
use polars::{datatypes::DataType, frame::DataFrame, lazy::frame::LazyFrame, prelude::*};
use rayon::prelude::*;
use std::{
    env,
    fmt::Write,
    sync::{Arc, Mutex},
    time::Duration,
    str::FromStr,
};

use crate::types::ContigMap;

pub fn create_subpileups(
    pileup: LazyFrame,
    contig_ids: Vec<String>,
    min_valid_read_coverage: u32,
    batches: usize,
) -> Result<Vec<DataFrame>> {
    info!("Creating subpileups");

    let tasks = contig_ids.len() as u64;
    let pb = ProgressBar::new(tasks);
    pb.set_draw_target(ProgressDrawTarget::stdout());
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len}",
        )
        .unwrap()
        .with_key("eta", |state: &ProgressState, w: &mut dyn Write| {
            write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap()
        })
        .progress_chars("#>-"),
    );
    pb.enable_steady_tick(Duration::from_millis(100));

    let mut subpileup_vec: Vec<DataFrame> = Vec::new();

    for contig_ids_batch in contig_ids.chunks(batches) {
        let pileup_batch = pileup
            .clone()
            .select([
                col("contig"),
                col("mod_type"),
                col("strand"),
                col("start"),
                col("N_valid_cov"),
                col("N_modified"),
            ])
            .filter(col("contig").is_in(lit(Series::new("contig".into(), contig_ids_batch))))
            .filter(col("N_valid_cov").gt_eq(lit(min_valid_read_coverage as i64)))
            .with_columns([(col("N_modified").cast(DataType::Float64)
                / col("N_valid_cov").cast(DataType::Float64))
            .alias("motif_mean")])
            .collect()
            .context("Error collecting pileup")?;

        if pileup_batch.is_empty() {
            error!("Pileup is empty after filtering");
            continue;
        }

        let mut subpileups = pileup_batch
            .partition_by(["contig"], true)
            .expect("Couldn't partition pileup");

        drop(pileup_batch);

        subpileup_vec.append(&mut subpileups);

        pb.inc(contig_ids_batch.len() as u64);
    }

    info!("Number of subpileups generated {:?}", subpileup_vec.len());
    Ok(subpileup_vec)
}

pub fn calculate_contig_read_methylation_pattern(
    contigs: ContigMap,
    subpileups: Vec<DataFrame>,
    motifs: Vec<Motif>,
    num_threads: usize,
) -> Result<DataFrame> {
    env::set_var("POLARS_MAX_THREADS", "1");

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Could not initialize threadpool");

    let tasks = subpileups.len() as u64;
    let pb = ProgressBar::new(tasks);
    pb.set_draw_target(ProgressDrawTarget::stdout());
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
    pb.enable_steady_tick(Duration::from_millis(100));

    let motifs = Arc::new(motifs);

    let schema = Schema::from_iter(vec![
        Field::new("contig".into(), DataType::String),
        Field::new("median".into(), DataType::Float64),
        Field::new("N_motif_obs".into(), DataType::UInt64),
        Field::new("mean_read_cov".into(), DataType::Float64),
        Field::new("motif".into(), DataType::String),
        Field::new("mod_type".into(), DataType::String),
        Field::new("mod_position".into(), DataType::Int32),
    ]);
    let empty_df = DataFrame::empty_with_schema(&schema);
    let results = Arc::new(Mutex::new(Vec::new()));

    subpileups.chunks(1000).for_each(|batch| {
        batch.par_iter().for_each(|subpileup| {
            let contig_id = match subpileup.column("contig") {
                Ok(column) => match column.get(0) {
                    Ok(value) => value.to_string().trim_matches('"').to_string(),
                    Err(e) => {
                        error!(
                            "Could not get first value of subpileup from contig column: {}",
                            e
                        );
                        return;
                    }
                },
                Err(e) => {
                    error!("Could not find contig column in subpileup: {}", e);
                    return;
                }
            };

            let contig_seq = match contigs.get(&contig_id) {
                Some(contig) => contig,
                None => {
                    error!("Could not find contig_id: {}", contig_id);
                    return;
                }
            };

            let mut local_read_methylation_df = empty_df.clone();

            for motif in motifs.iter() {
                let mod_type = motif.mod_type.to_pileup_code();

                let fwd_indices = find_motif_indices_in_contig(&contig_seq, &motif);

                let rev_indices =
                    find_motif_indices_in_contig(&contig_seq, &motif.reverse_complement());

                if fwd_indices.is_empty() && rev_indices.is_empty() {
                    continue;
                }

                let fwd_indices_series = Series::new("start".into(), fwd_indices);
                let rev_indices_series = Series::new("start".into(), rev_indices);

                let p_con = subpileup
                    .clone()
                    .lazy()
                    .filter(col("mod_type").eq(lit(mod_type)))
                    .filter(
                        (col("strand")
                            .eq(lit("+"))
                            .and(col("start").is_in(lit(fwd_indices_series.clone()))))
                        .or(col("strand")
                            .eq(lit("-"))
                            .and(col("start").is_in(lit(rev_indices_series.clone())))),
                    );

                let p_read_methylation = p_con
                    .group_by([col("contig")])
                    .agg([
                        (col("motif_mean").median()).alias("median"),
                        (col("contig").count()).alias("N_motif_obs"),
                        (col("N_valid_cov").cast(DataType::Float64).mean()).alias("mean_read_cov"),
                    ])
                    .with_columns([
                        (lit(motif.sequence.clone())).alias("motif"),
                        (lit(mod_type.to_string())).alias("mod_type"),
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

            let mut results_lock = results.lock().expect("Could not open results mutex");
            results_lock.push(local_read_methylation_df.lazy());
        });
        pb.inc(batch.len() as u64);
    });
    println!(""); //To stop log overwriting the progressbar

    let final_results: Vec<LazyFrame> = results.lock().unwrap().clone();

    let read_methylation_df = concat(&final_results, UnionArgs::default())
        .context("Unable to concatenate result dataframes.")?
        .collect()
        .context("Unable to collect result dataframes.")?;

    Ok(read_methylation_df)
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

            Motif::new(sequence, mod_type, mod_position).map_err(|e| anyhow!(e)).with_context(|| {
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
            calculate_contig_read_methylation_pattern(contig_map, subpileups, motifs, 1);

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
    fn test_create_subpileups() {
        let pileup = df!(
            "contig" => ["contig_3","contig_3","contig_3","contig_3","contig_3","contig_4","contig_4","contig_4","contig_4"],
            "strand" => ["+", "+", "+", "-", "-", "+", "+", "-", "-"],
            "mod_type" => ["a", "m", "a", "a", "a", "m", "a", "a", "a"],
            "start" => [6, 8, 12, 7, 13, 8, 12, 7, 13],
            "N_modified" => [20, 20, 5, 20 ,5, 2, 0, 20 ,5],
            "N_valid_cov" => [20 , 20, 20, 20, 20, 2, 2, 20, 20]
        )
        .expect("Could not init df")
        .lazy();

        let contig_ids = vec!["contig_3".to_string(), "contig_4".to_string()];

        let subpileups_1 = create_subpileups(pileup.clone(), contig_ids.clone(), 3 as u32, 2);

        assert_eq!(subpileups_1.len(), 2);

        for subpileup in subpileups_1 {
            let contig_id = match subpileup.column("contig") {
                Ok(column) => match column.get(0) {
                    Ok(value) => value.to_string().trim_matches('"').to_string(),
                    Err(e) => {
                        error!(
                            "Could not get first value of subpileup from contig column: {}",
                            e
                        );
                        return;
                    }
                },
                Err(e) => {
                    error!("Could not find contig column in subpileup: {}", e);
                    return;
                }
            };

            if contig_id == "contig_3" {
                assert_eq!(subpileup.shape().0, 5);
            } else {
                assert_eq!(subpileup.shape().0, 2);
            }
        }

        let subpileups_2 = create_subpileups(pileup, contig_ids, 1 as u32, 2);

        for subpileup in subpileups_2 {
            let contig_id = match subpileup.column("contig") {
                Ok(column) => match column.get(0) {
                    Ok(value) => value.to_string().trim_matches('"').to_string(),
                    Err(e) => {
                        error!(
                            "Could not get first value of subpileup from contig column: {}",
                            e
                        );
                        return;
                    }
                },
                Err(e) => {
                    error!("Could not find contig column in subpileup: {}", e);
                    return;
                }
            };

            if contig_id == "contig_3" {
                assert_eq!(subpileup.shape().0, 5);
            } else {
                assert_eq!(subpileup.shape().0, 4);
            }
        }
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
