use clap::Parser;
use motif::Motif;
use polars::{
    error::PolarsResult,
    lazy::frame::{LazyCsvReader, LazyFileListReader, LazyFrame},
    prelude::*,
};
use seq_io::fastq::{Reader, Record};
use std::path::Path;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long, required = true)]
    pileup: String,

    #[arg(short, long, required = true)]
    assembly: String,

    #[arg(short, long, default_value_t = 1)]
    threads: u8,

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
        "N_fail",
        "N_nocall",
    ];

    let lf_pileup = LazyCsvReader::new(path)
        .with_has_header(false)
        .with_separator(b'\t')
        .finish()?
        .rename(&old_column_names, &new_column_names, true);

    Ok(lf_pileup)
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

fn main() -> PolarsResult<()> {
    let args = Args::parse();

    let lf_pileup = load_pileup_lazy(&args.pileup)?;

    let mut fasta_reader = Reader::from_path(&args.assembly).unwrap();

    let motifs = match args.motifs {
        Some(motifs) => {
            println!("Motifs loaded");
            motifs
        }
        _ => panic!("No motifs found"),
    };

    let motifs = create_motifs(motifs);

    // let pileup = lf_pileup.collect()?;

    Ok(())
}