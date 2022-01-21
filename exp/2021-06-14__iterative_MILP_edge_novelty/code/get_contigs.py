__author__ = 'amane'

from sys import argv
import os
import time

import plasmids_preprocessing

sample_dir = '/home/aniket/python_scripts/Plasmids/data/unicycler_pipeline/'
contigs_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2019-05-07__iterative_MILP_unweighted/contigs/'

sample_id = argv[1]
assembly_file = sample_dir + 'sample_' + sample_id + '/assembly.gfa'
mapping_file = sample_dir + 'sample_' + sample_id + '/filtered_genes_to_contigs.csv'
seeds_file = sample_dir + 'sample_' + sample_id + '/seed_contigs.csv'

contigs_folder = contigs_dir + 'sample_' + sample_id
if not os.path.exists(contigs_folder):
    os.makedirs(contigs_folder)

contigs_file = open(os.path.join(contigs_folder,'putative_plasmids.fasta'),"w")

contigs_dict = {}
links_list = []

contigs_dict, links_list = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list)
contigs_dict = plasmids_preprocessing.get_gene_coverage(mapping_file, contigs_dict)

count = 0
for c in contigs_dict:
	count += 1
	print(count,c)
	contigs_file.write('>'+str(c)+'\t'+'Read_depth='+str(contigs_dict[c]['Read_depth'])+'\t'+'Gene_coverage='+str(contigs_dict[c]['Gene_coverage'])\
		+'\t'+'GC_cont='+str(contigs_dict[c]['GC_cont'])+'\t'+'Length='+str(contigs_dict[c]['Length'])+'\n')
	contigs_file.write(contigs_dict[c]['Sequence']+'\n')
	#print(c, contigs_dict[c]['Gene_coverage'])


		