use methylome::{ModType, Strand};
use serde::de::IgnoredAny;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct PileupRecord {
    // #[serde(rename = "0")]
    pub contig: String,

    // #[serde(rename = "1")]
    pub position: usize,

    // #[serde(rename = "2")]
    // #[serde(skip_deserializing)]
    _end: IgnoredAny,

    // #[serde(rename = "3")]
    pub mod_type: ModType,

    // #[serde(skip_deserializing)]
    // #[serde(rename = "4")]
    _score: IgnoredAny,

    // #[serde(rename = "5")]
    pub strand: Strand,

    // #[serde(skip_deserializing)]
    // #[serde(rename = "6")]
    _start2: IgnoredAny,

    // #[serde(skip_deserializing)]
    // #[serde(rename = "7")]
    _end2: IgnoredAny,

    // #[serde(skip_deserializing)]
    // #[serde(rename = "8")]
    _color: IgnoredAny,

    // #[serde(rename = "9")]
    pub n_valid_cov: u32,

    // #[serde(skip_deserializing)]
    // #[serde(rename = "10")]
    _percent_modified: IgnoredAny,

    // #[serde(rename = "11")]
    pub n_modified: u32,

    #[serde(skip)]
    // #[serde(rename = "12")]
    _n_canonical: IgnoredAny,

    #[serde(skip)]
    // #[serde(rename = "13")]
    _n_other_mod: IgnoredAny,

    #[serde(skip)]
    // #[serde(rename = "14")]
    _n_delete: IgnoredAny,

    #[serde(skip)]
    // #[serde(rename = "15")]
    _n_fail: IgnoredAny,

    #[serde(skip)]
    // #[serde(rename = "16")]
    _n_diff: IgnoredAny,

    #[serde(skip)]
    // #[serde(rename = "17")]
    _n_nocall: IgnoredAny,
}

#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::{BufReader, Write},
    };

    use csv::ReaderBuilder;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_deserialization() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);
        let mut rdr = ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(reader);

        for row in rdr.deserialize::<PileupRecord>() {
            let pileup_rec = row?;
            println!("{:?}", pileup_rec);
        }

        Ok(())
    }
}
