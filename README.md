# Methylation utils

Fast cli tool for extracting the read methylation of a pileup file. Supply the assembly, the pileup and the motifs of interest. The tool will:
 - Find motif occurences
 - Find the number of reads and mean read methylation at each position
 - calculate the median of mean methylated positions.

The return is a dataframe with:
 - contig
 - motif
 - mod_type
 - mod_position
 - median 
 - N_motif_obs
 - mean_read_cov


## Usage:
```bash
Cli tool for fast lookup in pileup for motif methylation

Usage: methylation_utils [OPTIONS] --pileup <PILEUP> --assembly <ASSEMBLY> --output <OUTPUT> --motifs <MOTIFS>...

Options:
  -p, --pileup <PILEUP>
          Path to pileup.
  -a, --assembly <ASSEMBLY>
          Path to assembly.
  -o, --output <OUTPUT>
          Path to output file. Must be .tsv.
  -t, --threads <THREADS>
          Number of parallel tasks. [default: 1]
  -m, --motifs <MOTIFS>...
          Supply chain of motifs as <motif>_<mod_type>_<mod_position>. Example: '-m GATC_a_1 RGATCY_a_2'
      --min-valid-read-coverage <MIN_VALID_READ_COVERAGE>
          Minimum valid read coverage for calculating methylation. [default: 3]
      --batches <BATCHES>
          Number of contigs to load at a time. Lower number will take longer but use less memory. 0 means no batching [default: 0]
  -h, --help
          Print help
  -V, --version
          Print version

```
