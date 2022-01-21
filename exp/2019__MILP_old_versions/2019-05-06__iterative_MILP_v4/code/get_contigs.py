__author__ = 'amane'

from sys import argv
import os
import time

sample_dir = '/home/aniket/python_scripts/Plasmids/data/unicycler_pipeline/'
contig_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2019-05-06__iterative_MILP_v4/contigs/'

sample_id = argv[1]
assembly_file = sample_dir + 'sample_' + sample_id + '/assembly.gfa'
mapping_file = sample_dir + 'sample_' + sample_id + '/filtered_genes_to_contigs.csv'
seeds_file = sample_dir + 'sample_' + sample_id + '/seed_contigs.csv'

contig_folder = contig_dir + 'sample_' + sample_id
if not os.path.exists(contig_folder):
    os.makedirs(contig_folder)

contigs_file = open(os.path.join(contigs_folder,'putative_plasmids.fasta'),"w")

string_list = read_file(assembly_file)
for line in string_list:
	line = line.split("\t")
	if line[0] == 'S':
		