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
2. A file mapping genes from a plasmid marker database to contigs (.tsv format)<br/>
Each line of this file contains the ids of the two sequences being mapped. In this case, we map genes to contigs, hence the file contains the gene name and the contig id. This is followed by the starting and ending positions of the contig sequence, denoting the region onto which the gene has been mapped.<br/>
Example:#Gene_name	Contig_ID	Start_Posn	End_Posn<br/>
AM349_24745	22	39217	38796<br/>
EC1303_109p00370	1	41519	41376<br/>
3. A file with a list of seed contigs (one entry per line)<br/>
4. Weight for the gene density term and GC content term in the objective function.

### Usage
```
python plasmids_iterative.py --ag assembly.gfa --map mapping.csv --seeds seed_contigs.csv \
				--out output_dir --alpha1 alpha_1 --alpha2 alpha_2 
```

## Output
### Contig chains
Lists of contigs and their orientation as they appear in the linear contig chain of all plasmid bins (both putative and questionable).<br/>
File name: contig_chains.csv, questionable_contig_chains.csv<br/>
Format: plasmid id;comma-separated list of oriented contigs<br/>
Example: <br/>
plasmid_0;23_0+,25_1-,10_0+<br/>
plasmid_0;23_1+,25_0<br/>


### Unique components
List of unique contig chains that are not subsets of other contig chains. </br>
File name: components.csv<br/>
Format: plasmid id;comma-separated list of contigs<br/>
Example: plasmid_0;23,25,10<br/>


### Sequences of contigs for unique components<br/>
Concatenated sequences of contigs in putative plasmid bins in FASTA format. The plasmid identifier is also used as the identifier of the FASTA entry. The additional information on each plasmid are provided in the deflines (e.g. seed_contig and gene_density.)<br/>
File names: putative_plasmids.fasta, questionable_plasmids.fasta<br/>
Format: FASTA with additional information in defline (as tab-separated list of <property>=<value> pairs)<br/>
