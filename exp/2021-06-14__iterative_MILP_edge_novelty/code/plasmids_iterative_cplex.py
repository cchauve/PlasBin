from __future__ import division
__author__ = 'amane'

#-------------------
#USAGE: 
#time python2.7 plasmids_iterative.py sample_id alpha_1 alpha_2 alpha_3 

from sys import argv
import collections
from random import *
import os
import time
import math

from docplex.mp.model import Model
from docplex.util.environment import get_environment
from docplex.mp.conflict_refiner import ConflictRefiner

import plasmids_prep_cplex

sample_dir = '/home/aniket/python_scripts/Plasmids/data/unicycler_pipeline/'
output_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2021-01-29__iterative_MILP_modified_graph/output/' 

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

#-----------------------------------------------
#Input files (Assembly graph, gene-contig mapping, seed contigs)
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
output_folder = output_dir + 'sample_' + sample_id + '/' + ratios + '/nplasmids_'+str(nplasmids) + '/MILP'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
#output_filename = argv[6]
output_filename = 'putative_plasmids.fasta'
questionable_filename = 'questionable_plasmids.fasta'
score_filename = 'MILP_objective.csv'
contigs_filename = 'contig_chains.csv'
q_contigs_filename = 'questionable_contig_chains.csv'

#-----------------------------------------------
#Main program
contigs_dict = {}
links_list = []
seeds_set = set()
rd_dict = {}
contigs_dict, links_list, rd_dict = plasmids_prep_cplex.get_data(assembly_file, contigs_dict, links_list, rd_dict)
seeds_set = plasmids_prep_cplex.get_seeds(seeds_file, seeds_set, rd_dict)




#print(seeds_set)

#GC_total = 0
ln_total = 0
n_contigs = 0
rd_graph = 0
for c in contigs_dict:
	if c in seeds_set:
		contigs_dict[c]['Seed'] = 1
	else:
		contigs_dict[c]['Seed'] = 0

	rd_graph += contigs_dict[c]['Read_depth']*contigs_dict[c]['Length']
	n_contigs += 1
	#GC_total += contigs_dict[c]['GC_cont']*contigs_dict[c]['Length']
	ln_total += contigs_dict[c]['Length']
rd_graph = 1
#print(rd_graph)	
#GC_mean = GC_total/ln_total		

#For consistency, a contig extremity should be a part of exactly one link for a specific plasmid.
#Here, for each extremity, we make a list of links that involve the extremity.
extr_dict = {}
for c in contigs_dict:
	ext1, ext2 = (c, 'h'), (c, 't')
	extr_dict[ext1] = []
	extr_dict[ext2] = []

	for link in links_list:
		if link[0] == ext1 or link[1] == ext1:
			extr_dict[ext1].append(link) 
		if link[0] == ext2 or link[1] == ext2:
			extr_dict[ext2].append(link)

contigs_dict = plasmids_prep_cplex.get_gene_coverage(mapping_file, contigs_dict, rd_dict)
#contigs_dict['1']['Gene_coverage'] = 0.8333
#contigs_dict['2']['Gene_coverage'] = 0.4163
#contigs_dict['3']['Gene_coverage'] = 0.9163
#contigs_dict['4']['Gene_coverage'] = 0.3333

UBD_rd = 0
UBD_GC = 0
UBD_gd = 0
for c in contigs_dict:
	#ln_total += contigs_dict[c]['Length']
	UBD_rd = max(UBD_rd, contigs_dict[c]['Read_depth'])
	UBD_GC = max(UBD_GC, contigs_dict[c]['GC_cont'])
	UBD_gd = max(UBD_gd, contigs_dict[c]['Gene_coverage'])
	print(c, contigs_dict[c]['Gene_coverage'], contigs_dict[c]['Read_depth'], contigs_dict[c]['Length'])	
print(UBD_rd, UBD_gd, UBD_GC, ln_total)


details_filename = 'details.csv'
details_file = open(os.path.join(output_folder, details_filename), "w")
details_file.write("Contig"+"\t"+'Gene_coverage'+ "\t"+'Read_depth'+ "\t"+'GC_cont'+"\t"+ "Length"+"\n")
for c in contigs_dict:
	details_file.write(c+"\t"+ str(contigs_dict[c]['Gene_coverage'])+ "\t"+str(contigs_dict[c]['Read_depth'])+ "\t"+str(contigs_dict[c]['GC_cont'])+"\t"+ str(contigs_dict[c]['Length'])+"\n")
	details_file.write("Link list:"+"\t")
	for link in links_list:
		if c == link[0][0] or c == link[1][0]:
			details_file.write(str(link)+"\t")
	details_file.write("\n")
			 
output_file = open(os.path.join(output_folder, output_filename), "w")
questionable_file = open(os.path.join(output_folder, questionable_filename), "w")
score_file = open(os.path.join(output_folder, score_filename), "w")
contigs_file = open(os.path.join(output_folder, contigs_filename), "w")
q_contigs_file = open(os.path.join(output_folder, q_contigs_filename), "w")
logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"w")

n_iter = 0
q_iter = 0

print(extr_dict)
print(links_list)

while len(seeds_set) > 0:
	print("\n\n\n\n\n")
	#-----------------------------------------------
	#Initializing the ILP
	m = Model("Plasmids")
	m.params.LogFile= os.path.join(output_folder,'m.log')
	m.setParam(GRB.Param.TimeLimit, 2400.0)

	#Minimum seed read depth check
	#-----------------------------------------------
	for x in seeds_set:
		print(x, contigs_dict[x]['Read_depth'])