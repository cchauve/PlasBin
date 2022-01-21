#!/bin/bash

for k in 1 2 3
do
	echo $k
	for ratio in '1.1.1' '1.1.0'
	do
		echo $ratio
		for id in 11 23 39 56 66 103 109 116 117 125	
		do
			echo $id					
			time python analyse_sample.py $id /home/amane/projects/rrg-chauvec/wg-anoph/Plasmids-Assembly/exp/2019-02-19__plasmid_MILP/output/sample_$id/$ratio/nplasmids_$k results.csv
		done	
	done
done
