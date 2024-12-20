use anyhow::Result;
use csv::StringRecord;
use methylome::{ModType, Strand};

use crate::data::{methylation::MethylationCoverage, MethylationRecord};

pub fn parse_to_methylation_record(
    contig: String,
    n_valid_cov: u32,
    record: &StringRecord,
) -> Result<MethylationRecord> {
    let position_str = record.get(1).expect("Missing position field");
    let position: usize = position_str.parse().expect("Error parsing position.");

    let mod_type_str = record.get(3).expect("Missing mod_type field.");
    let mod_type = ModType::from_str(mod_type_str)?;

    let strand_str = record.get(5).expect("Missing strand field.");
    let strand = Strand::from_str(strand_str)?;

    let n_modified_str = record.get(11).expect("Missing n_modified field.");
    let n_modified: u32 = n_modified_str
        .parse()
        .expect("Error parsing n_modified field");

    let methylation = MethylationCoverage::new(n_modified, n_valid_cov)?;

    let methylation_record =
        MethylationRecord::new(contig, position, strand, mod_type, methylation);

    Ok(methylation_record)
}
