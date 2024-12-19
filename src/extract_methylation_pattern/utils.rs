use anyhow::{anyhow, Context, Result};
use csv::StringRecord;
use methylome::{ModType, Strand};

use crate::data::{methylation::MethylationCoverage, MethylationRecord};

pub fn parse_to_methylation_record(
    contig: String,
    n_valid_cov: u32,
    record: &StringRecord,
    line_num: usize,
) -> Result<MethylationRecord> {
    let position_str = record
        .get(1)
        .ok_or_else(|| anyhow!("Missing position at line {}", line_num))?;
    let position: usize = position_str
        .parse()
        .with_context(|| format!("Invalid position at line {}", line_num))?;

    let mod_type_str = record
        .get(3)
        .ok_or_else(|| anyhow!("Missing mod_type at line {}", line_num))?;
    let mod_type = ModType::from_str(mod_type_str)?;

    let strand_str = record
        .get(5)
        .ok_or_else(|| anyhow!("Missing strand at line {}", line_num))?;
    let strand = Strand::from_str(strand_str)?;

    let n_modified_str = record
        .get(11)
        .ok_or_else(|| anyhow!("Missing n_modified at line {}", line_num))?;
    let n_modified: u32 = n_modified_str
        .parse()
        .with_context(|| format!("Invalid n_modified at line {}", line_num))?;

    let methylation = MethylationCoverage::new(n_modified, n_valid_cov)?;

    let methylation_record =
        MethylationRecord::new(contig, position, strand, mod_type, methylation);

    Ok(methylation_record)
}
