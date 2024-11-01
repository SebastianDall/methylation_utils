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


## Usage:
```bash
Usage: methylation_utils [OPTIONS] --pileup <PILEUP> --assembly <ASSEMBLY> --output <OUTPUT> --motifs <MOTIFS>...

Options:
  -p, --pileup <PILEUP>                                    
  -a, --assembly <ASSEMBLY>                                
  -o, --output <OUTPUT>                                    
  -t, --threads <THREADS>                                  [default: 1]
  -m, --motifs <MOTIFS>...                                 
      --min-valid-read-coverage <MIN_VALID_READ_COVERAGE>  [default: 3]
  -h, --help                                               Print help
  -V, --version                                            Print version

```
