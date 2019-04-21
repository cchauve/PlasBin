from __future__ import division
from gurobipy import *


def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

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
def update_contigs_dict(contigs_dict, line):
	c = get_id(line)
	seq = get_nucleotide_seq(line) 
	GC_cont = compute_GCratio(seq)
	ln = get_length(line)
	rd = get_read_depth(line)

	contigs_dict[c] = {}
	contigs_dict[c]['Sequence'] = seq
	contigs_dict[c]['Length'] = ln
	contigs_dict[c]['Read_depth'] = rd
	contigs_dict[c]['Seed'] = 0							#Default
	contigs_dict[c]['GC_cont'] = GC_cont
	contigs_dict[c]['Gene_coverage_intervals'] = []		#Default
	contigs_dict[c]['Gene_coverage'] = 0				#Default
	return contigs_dict

#A link is of the type: ((l1, e1),(l2, e2)) 
#where l1, l2 are adjacent links and e1, e2 are the connected link extremities
def get_link(line):
	c1, o1, c2, o2 = line[1], line[2], line[3], line[4]
	if o1 == '+':
		ext1 = 'h'
	else:
		ext1 = 't'
	if o2 == '+':
		ext2 = 't'
	else:
		ext2 = 'h'
	e = ((c1, ext1),(c2, ext2))	
	return e 

#Reads the assembly file line by line and forwards a line 
#to update_contigs_dict or get_link depending on the entry
def get_data(assembly_file, contigs_dict, links_list):
	string_list = read_file(assembly_file)
	count_s = 0
	count_l = 0
	for line in string_list:
		line = line.split("\t")
		if line[0] == 'S':
			contigs_dict = update_contigs_dict(contigs_dict, line)
		elif line[0] == 'L':
			e = get_link(line)
			links_list.append(e)
	return contigs_dict, links_list

#Reads the seed file and makes a list of seeds
def get_seeds(seeds_file, seeds_set):
	string_list = read_file(seeds_file)
	for line in string_list:
		line = line.split("\t")
		seeds_set.add(line[0])
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
def get_gene_coverage(mapping_file, contigs_dict):
	string_list = read_file(mapping_file)
	possible_seeds = []
	for line in string_list:
		line = line.split("\t")	
		qseqid, sseqid = line[0], line[1]
		sstart, send = line[8], line[9]
		if sseqid not in possible_seeds:
			possible_seeds.append(sseqid)
		if sseqid not in contigs_dict:
			print(sseqid, "not in contigs_dict")
		else:
			if int(sstart) > int(send):
				contigs_dict[sseqid]['Gene_coverage_intervals'].append((int(send), int(sstart)))
			else:
				contigs_dict[sseqid]['Gene_coverage_intervals'].append((int(sstart), int(send)))

	for sseqid in contigs_dict:
		union = get_union(contigs_dict[sseqid]['Gene_coverage_intervals'])
		ln = contigs_dict[sseqid]['Length']
		covered = 0
		for interval in union:
			covered += interval[1] - interval[0] + 1
		contigs_dict[sseqid]['Gene_coverage'] = covered/ln

	return contigs_dict	






#-----------------------------------------------------------	
#contigs[p][c] == 1 if contig 'c' belongs to plasmid 'p', else 0
def contig_vars(m, contigs_dict, contigs, contigs_ext, nplasmids):
	for p in range(nplasmids):
		contigs[p] = {}
		contigs_ext[p] = {}
		for c in contigs_dict:
			contigs[p][c] = m.addVar(vtype=GRB.BINARY, name='contig-'+c+'_plasmid-'+str(p))
			contigs_ext[p][(c, 'h')] = m.addVar(vtype=GRB.BINARY, name='contigext-'+c+'h'+'_plasmid-'+str(p))
			contigs_ext[p][(c, 't')] = m.addVar(vtype=GRB.BINARY, name='contigext-'+c+'t'+'_plasmid-'+str(p))
	return contigs, contigs_ext

#links[p][e] == 1 if link 'e' belongs to plasmid 'p', else 0
def link_vars(m, links_list, links, nplasmids):
	for p in range(nplasmids):
		links[p] = {}
		for e in links_list:		
			links[p][e] = m.addVar(vtype=GRB.BINARY, name='link-'+e[0][0]+e[0][1]+e[1][0]+e[1][1]+'_plasmid-'+str(p))
	return links

#-------------------------------------------------------------
def rd_vars(m, contigs_dict, mean_rd, rd_diff, counted_rd_diff, wtd_rd_diff, counted_wtd_rd_diff, nplasmids):
	for p in range(nplasmids):
		mean_rd[p] = {}
		rd_diff[p] = {}
		counted_rd_diff[p] = {}
		wtd_rd_diff[p] = {}
		counted_wtd_rd_diff[p] ={}
		mean_rd[p] = m.addVar(vtype=GRB.CONTINUOUS, name='mean-rd-plasmid-'+str(p))
		wtd_rd_diff[p] = m.addVar(vtype=GRB.CONTINUOUS, name='wtd-rd-diff-plasmid-'+str(p))
		for c in contigs_dict:	
			rd_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='rd-diff_contig-'+c+'_plasmid-'+str(p))
			counted_rd_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-rd-diff_contig-'+c+'_plasmid-'+str(p))
			counted_wtd_rd_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-wtd-rd-diff_contig-'+c+'_plasmid-'+str(p))
	return mean_rd, rd_diff, counted_rd_diff, wtd_rd_diff, counted_wtd_rd_diff

def GC_vars(m, contigs_dict, mean_GC, GC_diff, counted_GC_diff, wtd_GC_diff, counted_wtd_GC_diff, nplasmids):
	for p in range(nplasmids):
		mean_GC[p] = {}
		GC_diff[p] = {}
		counted_GC_diff[p] = {}
		wtd_GC_diff[p] = {}
		counted_wtd_GC_diff[p] ={}
		mean_GC[p] = m.addVar(vtype=GRB.CONTINUOUS, name='mean-GC-plasmid-'+str(p))
		wtd_GC_diff[p] = m.addVar(vtype=GRB.CONTINUOUS, name='wtd-GC-diff-plasmid-'+str(p))
		for c in contigs_dict:	
			GC_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='GC-diff_contig-'+c+'_plasmid-'+str(p))
			counted_GC_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-GC-diff_contig-'+c+'_plasmid-'+str(p))
			counted_wtd_GC_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-wtd-GC-diff_contig-'+c+'_plasmid-'+str(p))
	return mean_GC, GC_diff, counted_GC_diff, wtd_GC_diff, counted_wtd_GC_diff

def ln_vars(m, contigs_dict, counted_ln, nplasmids):
	for p in range(nplasmids):
		counted_ln[p] = {}	
		for c in contigs_dict:
			counted_ln[p][c] = m.addVar(vtype=GRB.INTEGER, name='counted-ln_contig-'+c+'_plasmid-'+str(p))
	return counted_ln		

def seed_vars(m, contigs_dict, counted_seed, nplasmids):
	for p in range(nplasmids):
		counted_seed[p] = {}
		for c in contigs_dict:		
			counted_seed[p][c] = m.addVar(vtype=GRB.INTEGER, name='counted-seed_contig-'+c+'_plasmid-'+str(p))
	return counted_seed

def gd_vars(m, contigs_dict, wtd_gd, counted_wtd_gd, nplasmids):
	for p in range(nplasmids):
		wtd_gd[p] = {}
		counted_wtd_gd[p] = {}
		wtd_gd[p] = m.addVar(vtype=GRB.INTEGER, name='wtd-gd_plasmid-'+str(p))
		for c in contigs_dict:			
			counted_wtd_gd[p][c] = m.addVar(vtype=GRB.INTEGER, name='counted-wtd-gd_contig-'+c+'_plasmid-'+str(p))
	return wtd_gd, counted_wtd_gd

#-----------------------------------------------------------------------		






			


