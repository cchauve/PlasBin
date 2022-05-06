from __future__ import division
from sys import argv
import os
import time
from random import randint
import math
import pandas as pd
import argparse

#-------------------

#USAGE: 
#time python generate_seeds.py --ag assembly.gfa --map mapping.csv --out output_dir \
#				 --gd_ratio gd_ratio --rd_ratio rd_ratio --max_len max_len

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
	rd = int(math.ceil(get_read_depth(line)))

	
	contigs_dict[c] = {}
	contigs_dict[c]['Sequence'] = seq
	contigs_dict[c]['Length'] = ln
	contigs_dict[c]['Read_depth'] = rd
	contigs_dict[c]['GC_cont'] = GC_cont
	contigs_dict[c]['Gene_coverage_intervals'] = []		#Default
	contigs_dict[c]['Gene_coverage'] = 0				#Default
	contigs_dict[c]['Density'] = 0				#Default

	return contigs_dict

#Reads the assembly file line by line and forwards a line 
#to update_contigs_dict or get_link depending on the entry
def get_data(assembly_file, contigs_dict):
	string_list = read_file(assembly_file)
	for line in string_list:
		line = line.split("\t")
		if line[0] == 'S':
			contigs_dict = update_contigs_dict(contigs_dict, line)
	return contigs_dict

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
	for line in string_list:
		line = line.split("\t")	
		qseqid, sseqid = line[0], line[1]
		sstart, send = line[8], line[9]

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

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--ag", help="Path to assembly graph file")
	parser.add_argument("--map", help="Path to gene to contig mapping file")
	parser.add_argument("--out", help="Path to output dir")
	parser.add_argument("--gd_ratio", nargs='?', const = 1, type=float, default = 1.5, help="Ratio between minimum seed gd and mean gd of assembly graph")
	parser.add_argument("--rd_ratio", nargs='?', const = 1, type=float, default = 0.3, help="Ratio between minimum seed read depth and median read depth of assembly graph")
	parser.add_argument("--max_len", nargs='?', const = 1, type=int, default = 1750000, help="Maximum length of seed contig")
	args = parser.parse_args()

	input_dir = args.input
	assembly_file = args.ag
	mapping_file = args.map

	contigs_dict = {}
	contigs_dict = get_data(assembly_file, contigs_dict)
	contigs_dict = get_gene_coverage(mapping_file, contigs_dict)

	contigs_df = pd.DataFrame.from_dict(contigs_dict).T
	mean_gd = contigs_df['Gene_coverage'].mean()
	med_rd = contigs_df['Read_depth'].median()

	min_seed_gd = float(args.gd_ratio)*mean_gd
	min_seed_rd = float(args.rd_ratio)*med_rd
	max_seed_len = int(args.max_len)

	seeds_file = open(os.path.join(input_dir, 'seed_contigs.csv'), "w")

	for c in contigs_dict:
		gd = contigs_dict[c]['Gene_coverage']
		rd = contigs_dict[c]['Read_depth']
		ln = contigs_dict[c]['Length']

		if gd >= min_seed_gd and rd >= min_seed_rd and ln <= max_seed_len:
			print(c)
			seeds_file.write(c+"\n")


