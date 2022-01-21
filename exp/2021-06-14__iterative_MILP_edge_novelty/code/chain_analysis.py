author__ = 'amane'

#-------------------

#USAGE: 
#time python2.7 chain_analysis.py sample_id alpha_1 alpha_2 alpha_3 

from gurobipy import *
from sys import argv
import os
import time
from random import randint
import math

import plasmids_preprocessing
import plasmids_postprocessing
import rmcircular

sample_dir = '/home/aniket/python_scripts/Plasmids/data/unicycler_pipeline/'
#sample_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2019-05-07__iterative_MILP_unweighted/'
#output_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2019-05-07__iterative_MILP_unweighted/output/' 


def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

#-----------------------------------------------
#Input files (Assembly graph, gene-contig mapping, seed contigs)
'''
assembly_file = argv[2]
mapping_file = argv[3]
seeds_file = argv[4]
'''

sample_id = argv[1]
assembly_file = sample_dir + 'sample_' + sample_id + '/assembly.gfa'
mapping_file = sample_dir + 'sample_' + sample_id + '/filtered_genes_to_contigs.csv'
seeds_file = sample_dir + 'sample_' + sample_id + '/seed_contigs.csv'

alpha1 = argv[2]
alpha2 = argv[3]
alpha3 = argv[4]

#-----------------------------------------------
#Output file 
#output_folder = argv[5]
ratios = alpha1 + '.' + alpha2 + '.' + alpha3
nplasmids = 1


#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#PATH DIFFERENCE ANALYSIS
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
details_file = '../output/sample_'+sample_id+'/' + ratios + '/nplasmids_'+str(nplasmids) + '/MILP/details.csv'
greedy_plasmids_file = '../contigs/sample_'+sample_id+'/greedy_contig_chains.csv'
MILP_plasmids_file = '../output/sample_'+sample_id+'/' + ratios + '/nplasmids_'+str(nplasmids) + '/MILP/contig_chains.csv'

details = read_file(details_file)
greedy_results = read_file(greedy_plasmids_file)
MILP_results= read_file(MILP_plasmids_file)


contigs_dict = {}
links_list = []

contigs_dict, links_list = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list)
contigs_dict = plasmids_preprocessing.get_gene_coverage(mapping_file, contigs_dict)

#print(contigs_dict)

for line in greedy_results:
	chain = line.split(';')[1]
	print(chain)
	chain = chain.split(',')
	for contig in chain:
		contig = contig[:-1]
		print(contig, contigs_dict[contig]['GC_cont'], contigs_dict[contig]['Gene_coverage'], contigs_dict[contig]['Length'])


