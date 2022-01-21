from __future__ import division
from gurobipy import *
import math


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
def update_contigs_dict(contigs_dict, line, rd_dict):
	c = get_id(line)
	seq = get_nucleotide_seq(line) 
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






#-----------------------------------------------------------	
#contigs[p][c] == 1 if contig 'c' belongs to plasmid 'p', else 0
def contig_vars(m, contigs_dict, contigs, degree, contigs_ext, ceil, nplasmids):
	for p in range(nplasmids):
		contigs[p] = {}
		degree[p] = {}
		contigs_ext[p] = {}
		ceil[p] = {}
		for c in contigs_dict:
			contigs[p][c] = m.addVar(vtype=GRB.BINARY, name='contig-'+c+'_plasmid-'+str(p))
			degree[p][c] = m.addVar(vtype=GRB.BINARY, name='degree-'+c+'_plasmid-'+str(p))
			contigs_ext[p][(c, 'h')] = m.addVar(vtype=GRB.BINARY, name='contigext-'+c+'h'+'_plasmid-'+str(p))
			contigs_ext[p][(c, 't')] = m.addVar(vtype=GRB.BINARY, name='contigext-'+c+'t'+'_plasmid-'+str(p))			
			ceil[p][c] = m.addVar(vtype=GRB.BINARY, name='ceil-'+c+'_plasmid-'+str(p))
	return contigs, degree, contigs_ext, ceil

#links[p][e] == 1 if link 'e' belongs to plasmid 'p', else 0
def link_vars(m, links_list, links, nplasmids):
	for p in range(nplasmids):
		links[p] = {}
		for e in links_list:		
			links[p][e] = m.addVar(vtype=GRB.BINARY, name='link-'+e[0][0]+e[0][1]+e[1][0]+e[1][1]+'_plasmid-'+str(p))
	return links

#-------------------------------------------------------------
def rd_vars(m, contigs_dict, rd_diff, counted_rd_diff, nplasmids):
	for p in range(nplasmids):
		rd_diff[p] = {}
		counted_rd_diff[p] = {}
		for c in contigs_dict:	
			rd_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='rd-diff_contig-'+c+'_plasmid-'+str(p))
			counted_rd_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-rd-diff_contig-'+c+'_plasmid-'+str(p))
	return rd_diff, counted_rd_diff

def GC_vars(m, contigs_dict, mean_GC, counted_GC_mean, GC_diff, counted_GC_diff, nplasmids):
	for p in range(nplasmids):
		mean_GC[p] = {}
		counted_GC_mean[p] = {}
		GC_diff[p] = {}
		counted_GC_diff[p] = {}
		mean_GC[p] = m.addVar(vtype=GRB.CONTINUOUS, name='mean-GC-plasmid-'+str(p))
		for c in contigs_dict:	
			GC_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='GC-diff_contig-'+c+'_plasmid-'+str(p))
			counted_GC_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-GC-diff_contig-'+c+'_plasmid-'+str(p))
			counted_GC_mean[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-GC-mean_contig-'+c+'_plasmid-'+str(p))
	return mean_GC, counted_GC_mean, GC_diff, counted_GC_diff

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
			counted_seed[p][c] = m.addVar(vtype=GRB.BINARY, name='counted-seed_contig-'+c+'_plasmid-'+str(p))
	return counted_seed

#def gd_vars(m, contigs_dict, wtd_gd, counted_wtd_gd, nplasmids):
#	for p in range(nplasmids):
#		wtd_gd[p] = {}
#		counted_wtd_gd[p] = {}
#		wtd_gd[p] = m.addVar(vtype=GRB.CONTINUOUS, name='wtd-gd_plasmid-'+str(p))
#		for c in contigs_dict:			
#			counted_wtd_gd[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-wtd-gd_contig-'+c+'_plasmid-'+str(p))
#	return wtd_gd, counted_wtd_gd

#-----------------------------------------------------------------------		






			


