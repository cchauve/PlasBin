#!/bin/bash

for k in 1
do
	echo $k
	for ratio in '1.1.1'
	do
		echo $ratio
		for id in 11 23 39 56 66 103 109 116 117 125	
		do
			echo $id					
			time python analyse_sample.py $id /home/aniket/python_scripts/Plasmids-Optimization/exp/2021-01-29__iterative_MILP_modified_graph/output/sample_$id/$ratio/nplasmids_$k results.csv
		done	
	done
done
