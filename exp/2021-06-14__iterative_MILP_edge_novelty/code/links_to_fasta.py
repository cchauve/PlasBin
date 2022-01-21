from __future__ import division

__author__ = 'amane'

from sys import argv
import os
import time
import math

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

#Obtain set of contigs in a plasmid
def get_contigs(plasmid):
	plasmid = plasmid.split(',')
	contigs = set()
	for x in plasmid:
		c = x.split('_')[0]
		contigs.add(c)
	return contigs

#Storing contig details
#-----------------------------------------------
#Stores the id of the contig
def get_id(line):
	return line[1]
#Stores the nucleotide sequence of the contig
def get_nucleotide_seq(line):
	#print(line[2])
	return line[2]		
#Computes GC ratio: counts no. of 'G'/'C' occurences in the sequence and divide by the sequence length.
def compute_GCratio(seq):
	GC = 0
	ln_seq = 0
	for nucl in seq:
		if nucl == 'G' or nucl == 'C':
			GC += 1
		ln_seq += 1
	return GC/ln_seq
#Stores the length of the sequence
def get_length(line):
	return int(line[3].split(':')[2])
#Stores the read depth of the contig
def get_read_depth(line):
	return float(line[4].split(':')[2])		

#Takes a contig from the assembly file and initiates an entry in the contigs_dict
#Each contig is tagged with the following attributes:
#1. Length of the contig (int)
#2. Overall read depth of the contig (float)
#3. Indication if the contig is a seed (binary)
#4. GC content of the contig (float)
#5. Gene coverage intervals (list of pairs)
#6. Gene coverage (float)
def update_contigs_dict(contigs_dict, line, rd_dict):
	c = get_id(line)
	seq = get_nucleotide_seq(line) 
	#print(seq)
	GC_cont = compute_GCratio(seq)
	ln = get_length(line)
	rd = int(math.ceil(get_read_depth(line)))

	for i in range(rd):
		contigs_dict[c+"_"+str(i)] = {}
		contigs_dict[c+"_"+str(i)]['Sequence'] = seq
		contigs_dict[c+"_"+str(i)]['Length'] = ln
		contigs_dict[c+"_"+str(i)]['Read_depth'] = 1
		contigs_dict[c+"_"+str(i)]['Seed'] = 0							#Default
		contigs_dict[c+"_"+str(i)]['GC_cont'] = GC_cont
		contigs_dict[c+"_"+str(i)]['Gene_coverage_intervals'] = []		#Default
		contigs_dict[c+"_"+str(i)]['Gene_coverage'] = 0				#Default
		contigs_dict[c+"_"+str(i)]['Density'] = 0				#Default
	rd_dict[c] = rd	
	return contigs_dict,rd_dict	

#A link is of the type: ((l1, e1),(l2, e2)) 
#where l1, l2 are adjacent links and e1, e2 are the connected link extremities
def get_link(line, rd_dict, links_list, edge_fam_dict):
	c1, o1, c2, o2 = line[1], line[2], line[3], line[4]
	if o1 == '+':
		ext1 = 'h'
	else:
		ext1 = 't'
	if o2 == '+':
		ext2 = 't'
	else:
		ext2 = 'h'
	rd1 = rd_dict[c1]
	rd2 = rd_dict[c2]

	edge = ((c1, ext1),(c2, ext2))
	edge_fam = []
	for i in range(rd1):
		for j in range(rd2):
			e = ((c1+"_"+str(i), ext1),(c2+"_"+str(j), ext2))
			links_list.append(e)
			edge_fam.append(e)	
	if edge not in edge_fam_dict:
		edge_fam_dict[edge] = edge_fam		
	return links_list, edge_fam_dict

#Reads the assembly file line by line and forwards a line 
#to update_contigs_dict or get_link depending on the entry
def get_data(assembly_file, contigs_dict, links_list,rd_dict, edge_fam_dict):
	string_list = read_file(assembly_file)
	count_s = 0
	count_l = 0
	repeated_contigs = set()
	for line in string_list:
		line = line.split("\t")
		if line[0] == 'S':
			contigs_dict, rd_dict = update_contigs_dict(contigs_dict, line, rd_dict)
		elif line[0] == 'L':
			links_list, edge_fam_dict = get_link(line, rd_dict, links_list, edge_fam_dict)
	return contigs_dict, links_list, rd_dict, edge_fam_dict	

#Reads the seed file and makes a list of seeds
def get_seeds(seeds_file, seeds_set, rd_dict):
	string_list = read_file(seeds_file)
	for line in string_list:
		line = line.split("\t")
		for i in range(rd_dict[line[0]]):
			seeds_set.add(line[0]+"_"+str(i))
	return seeds_set

#Takes the gene covering intervals for a contig and finds their union
#The length of the union is used to compute gene coverage
def get_union(intervals):
	union = []
	for begin, end in sorted(intervals):
		if union and union[-1][1] >= begin - 1:
			union[-1][1] = max(union[-1][1], end)
		else:
			union.append([begin, end])
	return union		

#Computes the gene coverage for each contig
def get_gene_coverage(mapping_file, contigs_dict, rd_dict):
	string_list = read_file(mapping_file)
	possible_seeds = []
	for line in string_list:
		line = line.split("\t")	
		qseqid, sseqid = line[0], line[1]
		sstart, send = line[8], line[9]

		if sseqid not in possible_seeds:
			possible_seeds.append(sseqid)
		if sseqid not in rd_dict:
			print(sseqid, "not in contigs_dict or rd_dict")
		else:
			rd = rd_dict[sseqid]
			for i in range(rd):
				new_sseqid = sseqid+'_'+str(i)
				if int(sstart) > int(send):
					contigs_dict[new_sseqid]['Gene_coverage_intervals'].append((int(send), int(sstart)))
				else:
					contigs_dict[new_sseqid]['Gene_coverage_intervals'].append((int(sstart), int(send)))

	for sseqid in contigs_dict:
		union = get_union(contigs_dict[sseqid]['Gene_coverage_intervals'])
		ln = contigs_dict[sseqid]['Length']
		covered = 0
		for interval in union:
			covered += interval[1] - interval[0] + 1
		contigs_dict[sseqid]['Gene_coverage'] = covered/ln
		if contigs_dict[sseqid]['Gene_coverage'] > 0:
			contigs_dict[sseqid]['Density'] = 1

	return contigs_dict			
			
sample_id = argv[1]

output_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2021-06-14__iterative_MILP_edge_novelty/output/'
sample_dir = '/home/aniket/python_scripts/Plasmids/data/unicycler_pipeline/'
assembly_file = sample_dir + 'sample_' + sample_id + '/assembly.gfa'
mapping_file = sample_dir + 'sample_' + sample_id + '/filtered_genes_to_contigs.csv'
seeds_file = sample_dir + 'sample_' + sample_id + '/seed_contigs.csv'

links_file = output_dir + 'sample_' + sample_id + '/1.1.1/nplasmids_1/MILP/links.csv'
details_file = output_dir + 'sample_' + sample_id + '/1.1.1/nplasmids_1/MILP/details.csv'

links_new = output_dir + 'sample_' + sample_id + '/1.1.1/nplasmids_1/MILP/final_chains.csv'
output_file = output_dir + 'sample_' + sample_id + '/1.1.1/nplasmids_1/MILP/putative_plasmids.fasta'

		
#Getting details
contigs_dict = {}
links_list = []
rd_dict = {}
edge_fam_dict = {}
seeds_set = set()
contigs_dict, links_list, rd_dict, edge_fam_dict = get_data(assembly_file, contigs_dict, links_list, rd_dict, edge_fam_dict)
seeds_set = get_seeds(seeds_file, seeds_set, rd_dict)
contigs_dict = get_gene_coverage(mapping_file, contigs_dict, rd_dict)

string_list = read_file(links_file)
plasmids_list = [x for x in string_list if x[-1] != ';'] #Individual components as separate plasmids

contigs_list = []
count = 0
for plasmid in plasmids_list: #For every individual plasmid component
	contigs = get_contigs(plasmid)

	temp_list = []
	flag = 0 #flag = 1 if set of contigs of current plasmid is subset of set of contigs of any previous plasmids
	for contig_set in contigs_list:
		if contig_set[0].issubset(contigs) == False: 
			temp_list.append(contig_set)
		elif contig_set[0] == contigs:
			temp_list.append(contig_set)
		if flag == 0 and contigs.issubset(contig_set[0]):
			flag = 1
	if flag == 0:
		temp_list.append((contigs,count))
	contigs_list = temp_list	
	count += 1			

# Returns that starting and ending point (index) of the sublist, if it exists, otherwise 'None'.
def find_sub_list(sub_list, in_list):
	print("Inside find sub list")
	print(sub_list[0])
	print(type(in_list))
	sub_list_length = len(sub_list[0])
	flag = 0
	for i in range(len(in_list)-sub_list_length):
		if sub_list[0] == in_list[i:i+sub_list_length]:
			flag = 1
			return (i,i+sub_list_length)
	if flag == 0:
		return None	       


# Removes the sublist, if it exists and returns a new list, otherwise returns the old list.
def remove_sub_list(sub_list, in_list):
	indices = find_sub_list(sub_list, in_list)
	print("F1", indices)
	if not indices is None:
		return in_list[0:indices[0]] + in_list[indices[1]:]
	else:
		return in_list	




for pair in contigs_list:
	idx = pair[1]
	temp_len = 0
	temp_seq = []
	gene_free_len = []
	gene_free_seq = []
	print("\n")
	print(idx, plasmids_list[idx])
	chain = plasmids_list[idx].split(',')
	plasmids_list[idx] = chain
	
	for contig in chain:
		gd = contigs_dict[contig[:-1]]['Gene_coverage']
		ln = contigs_dict[contig[:-1]]['Length']
		#print(contig[:-1], gd, ln)

		if gd >= 0.3:
			if temp_len > 11000:
				gene_free_len.append(temp_len)
				gene_free_seq.append(temp_seq)
			temp_len = 0
			temp_seq = []
		else:
			temp_len += ln
			temp_seq.append(contig)

	if len(gene_free_seq) == 0:
		print("No gene free stretches")		

	for x in range(len(gene_free_seq)):
		print(idx, gene_free_len[x], gene_free_seq[x])
		print(remove_sub_list(gene_free_seq,chain))
		plasmids_list[idx] = remove_sub_list(gene_free_seq,chain)
					







#Printing unique plasmids to output file
output_filename = open(os.path.join(output_file), "w")
links_new_filename = open(os.path.join(links_new), "w")
count = 1
for contig_set in contigs_list:
	c = contig_set[1]
	#print(contig_set)
	
	#plasmids_list[c] = plasmids_list[c].split(',')
	#print(plasmids_list[c])	
	chain = plasmids_list[c]
	#links_new.write(">plasmid_"+str(count)+"\t"+str(chain))
	print("\n")
	print(count, chain)
	if len(chain) > 1 and chain[0] == chain[-1]:
		chain = chain[:-1]
	length = 0
	gd = 0
	rd = 1
	seq = ''
	links_new_filename.write("plasmid_"+str(count)+";")
	for contig in chain:
		print(contig)
		length += contigs_dict[contig[:-1]]['Length']
		gd += contigs_dict[contig[:-1]]['Gene_coverage']*contigs_dict[contig[:-1]]['Length']
		seq += contigs_dict[contig[:-1]]['Sequence']
		links_new_filename.write(str(contig)+",")
	gd = gd/length
	if gd >= 0.3 and length >= 1500: 	
		links_new_filename.write("\n")

		output_filename.write(">plasmid_"+str(count)+"\t"+"length="+str(length)+"\t"+"gene_density="+str(gd)+"\t"+"mean_read_depth="+str(1)+"\n")
		output_filename.write(seq+"\n")
	else:
		print("QUESTIONABLE")
		print("plasmid_",count, length, gd)
		print("\n")

	
	

	count+=1

