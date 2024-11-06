use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[arg(short, long, required = true)]
    pub pileup: String,

    #[arg(short, long, required = true)]
    pub assembly: String,

    #[arg(short, long, required = true)]
    pub output: String,

    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,

    #[arg(short, long, required = true, num_args(1..))]
    pub motifs: Option<Vec<String>>,

    #[arg(long, default_value_t = 3)]
    pub min_valid_read_coverage: u32,
}
