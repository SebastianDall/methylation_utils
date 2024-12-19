use clap::Parser;

#[derive(Parser, Debug, Clone)]
pub struct MethylationPatternArgs {
    #[arg(short, long, required = true, help = "Path to pileup.")]
    pub pileup: String,

    #[arg(short, long, required = true, help = "Path to assembly.")]
    pub assembly: String,

    #[arg(
        short,
        long,
        required = true,
        help = "Path to output file. Must be .tsv."
    )]
    pub output: String,

    #[arg(short, long, default_value_t = 1, help = "Number of parallel tasks.")]
    pub threads: usize,

    #[arg(short, long, required = true, num_args(1..), help = "Supply chain of motifs as <motif>_<mod_type>_<mod_position>. Example: '-m GATC_a_1 RGATCY_a_2'")]
    pub motifs: Option<Vec<String>>,

    #[arg(
        long,
        default_value_t = 3,
        help = "Minimum valid read coverage for calculating methylation."
    )]
    pub min_valid_read_coverage: u32,

    #[arg(
        long,
        default_value_t = 10000,
        help = "Number of contigs to load at a time. Lower number will take longer but use less memory."
    )]
    pub batches: usize,
}
