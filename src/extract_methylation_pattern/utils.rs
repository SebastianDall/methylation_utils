use anyhow::{anyhow, Result};
use csv::StringRecord;
use methylome::{ModType, Strand};

use crate::data::{methylation::MethylationCoverage, MethylationRecord};

pub fn parse_to_methylation_record(
    contig: String,
    n_valid_cov: u32,
    record: &StringRecord,
) -> Result<MethylationRecord> {
    let position: usize = record
        .get(1)
        .ok_or_else(|| anyhow!("Missing position field."))?
        .parse()
        .map_err(|| anyhow!("Invalid position field"));

    let mod_type: ModType = record
        .get(3)
        .ok_or_else(|| anyhow!("Missing modification type field."))?
        .parse()?;

    let strand: Strand = record
        .get(5)
        .ok_or_else(|| anyhow!("Missing strand field"))?
        .parse()?;

    let n_modified: u32 = record
        .get(11)
        .ok_or_else(|| anyhow!("Missing n_modified field."))?
        .parse()
        .map_err(|| "Invalid n_modified field")?;

    let methylation = MethylationCoverage::new(n_modified, n_valid_cov)?;

    let methylation_record =
        MethylationRecord::new(contig, position, strand, mod_type, methylation);

    Ok(methylation_record)
}
