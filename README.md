# epimetheus
`epimetheus` is a CLI tool for converting a pileup file ([modkit](https://github.com/nanoporetech/modkit)), where reads with modification calls are mapped to an assembly, to motif methylation instead.

To use `epimetheus` it is therefore necessary to know which motifs to look for beforehand (epimetheus meaning hindsight or afterthought). [Nanomotif](https://github.com/MicrobialDarkMatter/nanomotif) can find motifs in metagenomic data using Nanopore methylation calls.



## Usage:
```bash
Cli tool for fast lookup in pileup for motif methylation

Usage: epimetheus <COMMAND>

Commands:
  methylation-pattern  
  help                 Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```

### methylation pattern

Efficient processing of a pileup file for finding the read methylation degree of a motif for all contigs. Supply the assembly, the pileup and the motifs of interest. The tool will:
 - Find motif occurences
 - Find the number of reads and mean read methylation at each position
 - calculate the median of mean methylated positions.

The return is a dataframe with:
- contig: The contig id
- motif: The motif sequence
- mod_type: The modification type [6mA, 5mC, 4mC (as pileup codes)]
- mod_position: The modification position in the motif sequence
- median: The median motif read methyalation
- mean_read_cov: The mean read coverage for positions used in the median calculation
- N_motif_obs: The number of motifs with methylation information above `min-valid-read-coverage`
- motif_occurences_total: The total of occurences of the motif sequence in the contig.

```bash
Usage: epimetheus methylation-pattern [OPTIONS] --pileup <PILEUP> --assembly <ASSEMBLY> --output <OUTPUT> --motifs <MOTIFS>...

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
          Number of contigs to process at a time. Higher number will use more RAM. [default: 3000]
  -h, --help
          Print help
```
