from sys import argv
import os

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

sample_id = argv[1]

output_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2021-01-29__iterative_MILP_modified_graph/output/'

output_file = output_dir + 'sample_' + sample_id + '/1.1.1/nplasmids_1/MILP/putative_plasmids.fasta'