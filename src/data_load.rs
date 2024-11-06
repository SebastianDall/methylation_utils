use polars::{
    error::PolarsResult,
    lazy::frame::{LazyCsvReader, LazyFileListReader, LazyFrame},
};
use seq_io::fasta::{Reader, Record};
use std::{path::Path};

use crate::types::ContigMap;

pub fn load_pileup_lazy<P: AsRef<Path>>(path: P) -> PolarsResult<LazyFrame> {
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

pub fn load_contigs<P: AsRef<Path>>(path: P) -> Result<ContigMap, Box<dyn std::error::Error>> {
    let mut fasta_reader = Reader::from_path(path).unwrap();
    let mut contigs = ContigMap::new();

    while let Some(record_result) = fasta_reader.next() {
        let record = record_result?;

        let id = record.id()?.to_string();
        let seq = String::from_utf8(record.owned_seq())?.to_string();

        contigs.insert(id, seq);
    }
    Ok(contigs)
}
