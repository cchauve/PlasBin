#!/bin/bash

for k in 1
do
	echo $k
	for ratio in '1.1.1'
	do
		echo $ratio
		for id in {62,63,64,65,66,76,85,86,87}	
		do
			echo $id		
			
			#mv /home/aniket/python_scripts/Plasmids-Optimization/exp/2021-06-14__iterative_MILP_edge_novelty/output/sample_$id/$ratio/nplasmids_$k/MILP/putative_plasmids.fasta /home/aniket/python_scripts/Plasmids-Optimization/exp/2021-06-14__iterative_MILP_edge_novelty/output/sample_$id/$ratio/nplasmids_$k/MILP/MILP_putative_plasmids.fasta
			#cp /home/aniket/python_scripts/Plasmid_Comparison/HyAsP_results/sample_$id/outputs/putative_plasmids.fasta /home/aniket/python_scripts/Plasmids-Optimization/exp/2021-06-14__iterative_MILP_edge_novelty/output/sample_$id/$ratio/nplasmids_$k/MILP/putative_plasmids.fasta
			time python analyse_sample.py $id /home/aniket/python_scripts/Plasmids-Optimization/exp/2021-06-14__iterative_MILP_edge_novelty/output/sample_$id/$ratio/nplasmids_$k results.csv
		done	
	done
done
