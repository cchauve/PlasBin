__author__ = 'amane'

from sys import argv
import os
import time

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list



contigs_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2019-05-07__iterative_MILP_unweighted/contigs/'
sample_id = argv[1]
putative_map = contigs_dir + 'sample_' +sample_id + '/putative_map.csv'

chain_file = open(os.path.join(contigs_dir, 'sample_'+sample_id,'contig_chains.csv'),"w")

string_list = read_file(putative_map)

contig_intervals = {}

for line in string_list:
	line = line.split("\t")
	contig, plasmid_id, coverage, sstart, send = line[0], line[1], float(line[2]), int(line[8]), int(line[9])
	if sstart < send:
		contig = contig+'+'
		interval = [sstart, send]
	else:
		contig = contig+'-'
		interval = [send, sstart]

	if plasmid_id not in contig_intervals:
		contig_intervals[plasmid_id] = []
		contig_intervals[plasmid_id].append((contig, interval))
	else:
		contig_intervals[plasmid_id].append((contig, interval))

possible_chains = {}

for plasmid in contig_intervals:	
	#print(plasmid)
	contig_intervals[plasmid] = sorted(contig_intervals[plasmid], key=lambda x: x[1][0])

	possible_chains[plasmid] = {}
	send = 0

	for pair in contig_intervals[plasmid]:
		#print("\n")
		#print(pair)
		#print(len(possible_chains[plasmid]))		
		sstart = pair[1][0]
		chain_found = 0

		if len(possible_chains[plasmid]) == 0:
			possible_chains[plasmid][0] = [pair]
			#print(len(possible_chains[plasmid]))
		else:
			for x in possible_chains[plasmid]:
				#print(x)
				#print(possible_chains)
				chain = possible_chains[plasmid][x]
				send = chain[-1][1][1]
				if sstart > send:
					chain_found = 1
					chain.append(pair)
					#print(chain)
					possible_chains[plasmid][x] = chain

			if chain_found == 0:
				n = len(possible_chains[plasmid])
				possible_chains[plasmid][n] = [pair]	

	#print("\n")
	#print(len(possible_chains[plasmid]))	

	opt_x = 0
	max_length = 0
	for x in possible_chains[plasmid]:
		interval_length = 0
		for pair in possible_chains[plasmid][x]:
			sstart, send = pair[1][0], pair[1][1]
			interval_length += (send - sstart + 1)
		if interval_length > max_length:
			max_length = interval_length
			opt_x = x

	chain_file.write(plasmid+';')
	for pair in possible_chains[plasmid][opt_x]:
		if pair == possible_chains[plasmid][opt_x][-1]:
			print(pair[0])
			chain_file.write(pair[0]+'\n')
		else:
			print(pair[0])
			chain_file.write(pair[0]+',')					

				

#print(contig_intervals)					

			


