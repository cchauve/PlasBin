# PlasBin
PlasBin is a tool that uses a mixed integer linear programming (MILP) based tool for plasmid binning. It uses a hybrid approach that incorporates ideas from both de novo and reference-based methods to identify plasmid bins. PlasBin considers several features associated with a contig such as
its %GC content, sequencing coverage, length and plasmid gene density.

## Overview
The code/ directory contains the source code of PlasBin. The requirements to run PlasBin have been listed below. The results/ directory contains the results of PlasBin, plasmidSPAdes, MOB-recon, HyAsP and gplas.

## Requirements
The following are required to run HyAsP
1. Python (Version 2.7+; packages: random, math, sys)
2. Gurobi solver (Version 9.1.2+)

## Input
1. A file containing the details of the assembly graph (.gfa format), 
2. A file mapping genes from a plasmid marker database to contigs (.csv format)
3. A file with a list of seed contigs (one entry per line) 
4. Weight for the gene density term and GC content term in the objective function.

### Usage
```
python2.7 plasbin_iterative.py assembly.gfa genes_to_contigs.csv seed_contigs.txt alpha_1 alpha_2 alpha_3
```

## Output
### Contig chains
Lists of contigs and their orientation as they appear in the linear contig chain of all plasmid bins (both putative and questionable).
File name: putative_contig_chains.csv, questionable_contig_chains.csv
Format: <plasmid id>;<comma-separated list of contigs with orientation>
Example: plasmid_0;23+,25-,10+

### Sequences of contig chains
Concatenated sequences of contigs in putative plasmid bins in FASTA format. The plasmid identifier is also used as the identifier of the FASTA entry. The additional information on each plasmid are provided in the deflines (e.g. seed_contig and gene_density.)
File names: putative_plasmids.fasta, questionable_plasmids.fasta
Format: FASTA with additional information in defline (as tab-separated list of <property>=<value> pairs)
