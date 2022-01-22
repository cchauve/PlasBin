# Plasmids-Optimization
This repository contains the details and experimental results of the Plasmid MILP tool. 
It also contains the results of five other methods designed for plasmid binning and assembly. 
These methods include HyAsP, MOB-recon, plasmidSPAdes, gplas and SCAPP.

## Data
We use the data set from Robertson and Nash, 2018 (previously used for the HyAsP experiments). 
The data set contains 133 bacterial samples. As some of the methods we use are reference-based, we split this data set into two parts - a reference set and a test set. 
We carry out our experiments on the test set consisting of 66 of the 133 samples.

*Question.* How was the plit done?   
*Comment.* If the reference set is used, we should outline how.

## Experiment details

## MILP
The MILP is a hybrid approach aimed at extracting cycles or paths from the assembly graph that represent plasmids *(fragments?)*. 
The objective function is inspired by the greedy heuristic in that it tries to maximize the gene density of plasmids and maintain a uniform coverage and GC content. The MILP iteratively outputs a plasmid that contains at least one seed contig.

*Comment.* Here we should make it clear if we bin or assemble. We need to provide a slightly more detailed overview of the MILP method.

### Input
1. Assembly graph, 2. Gene density of contigs, 3. Read depth of each contig, 4. GC content of each contig, 5. Seed eligbility of each contig

*Comment.* We will need to specify how the contig features are obtained.

### Output
1. Contig chains representing plasmids, 2. Fasta file containing the sequences of assembled *(assembled or binned?)* plasmids.

Older versions of the MILP can be found in the folder exp/2019__MILP_old_versions. 
The two latest versions are located in the exp folder with their own directories.
### 2021-01-29__iterative_MILP_modified_graph
This version executed the changes pertaining to the input assembly graph. 
These changes were motivated by the need to visit the same node multiple times, something that wasn't easily achievable in previous versions. 
As a result, a contig with coverage $k$ is replaced by $k$ copies of the contig. 
The edges incident on the original contig are now incident on all $k$ copies of the contig. 
Since every new contig copy has coverage 1, this also allows us to remove the read depth uniformity term in the MILP

*Comment.* This is hard to follow without a description of the initial MILP.

### 2021-06-14__iterative_MILP_edge_novelty
This version of the MILP included edge novelty constraints. 
These constraints allow the MILP to choose solutions that allow multiple copies of the same contig but not similar edges (If $(a_1,b_1)$ is chosen then no other edge $(a_x,b_y)$ will be chosen. 

## HyAsP
HyAsP is a hybrid approach that uses a greedy heuristic to extract cycles or paths from an assembly graph. 
It chooses a seed contig as a starting point and extends the path based on an objective function, similar to the MILP. 
However, in case of HyAsP, it can only extend based on the immediate neighbors of the current end contig.  
### Input
1. Assembly graph, 2. Gene density of contigs, 3. Read depth of each contig, 4. GC content of each contig, 5. Seed eligbility of each contig
### Output
1. Contig chains representing plasmids, 2. Fasta file containing the sequences of assembled plasmids.

## plasmidSPAdes


## MOB-recon (from MOB-suite)


## gplas


## SCAPP


## Documents
The analysis of the latest results is provided in the notebook results_analysis.ipynb, located in the folder exp/2021-06-14__iterative_MILP_edge_novelty/doc. 
The notebook contains the precision, recall and F1 scores for the results of the MILP, HyAsP, MOB-recon and plasmidSPAdes. 
It also contains a comparative analysis between the MILP and HyAsP.
