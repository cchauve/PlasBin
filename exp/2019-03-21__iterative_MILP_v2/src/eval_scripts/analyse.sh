#!/bin/bash

for k in 1
do
	echo $k
	for ratio in '1.1.1' '1.1.0' '1.2.1' '1.5.1'
	do
		echo $ratio
		for id in 11 23 39 56 66 103 109 116 117 125	
		do
			echo $id					
			time python analyse_sample.py $id /home/amane/projects/rrg-chauvec/wg-anoph/Plasmids-Optimization/exp/2019-03-21__iterative_MILP_v2/output/sample_$id/$ratio/nplasmids_$k results.csv
		done	
	done
done
