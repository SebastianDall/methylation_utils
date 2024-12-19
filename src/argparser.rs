use crate::extract_methylation_pattern::args::MethylationPatternArgs;
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    MethylationPattern(MethylationPatternArgs),
}
