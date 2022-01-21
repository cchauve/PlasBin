from sys import argv
import os
import time
from random import randint
import math

ILP_dir = "../output/"
greedy_dir = "../contigs/"

sample_id = argv[1]

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

details = ILP_dir + 'sample_' + sample_id + '/1.1.1/nplasmids_1/MILP/details.csv'
ILP_plasmids = ILP_dir + 'sample_' + sample_id + '/1.1.1/nplasmids_1/MILP/contig_chains.csv'
greedy_file = greedy_dir + 'sample_' + sample_id + '/greedy_contig_chains.csv'
real_file = greedy_dir + 'sample_' + sample_id + '/contig_chains.csv'

def eval_obj(seq, c_dict):
	total_len = 0
	total_score = 0
	avg_gcl = 0
	for c in seq:
		#c = c[:-1]
		gc = c_dict[c]['gc']
		cov = c_dict[c]['cov']
		length = c_dict[c]['length']

		avg_gcl += gc*length
		#score = (cov - gc)*length
		#total_score += score
		total_len += length 

	avg_gcl = avg_gcl/total_len
	
	for c in seq:
		gc = c_dict[c]['gc']
		cov = c_dict[c]['cov']
		length = c_dict[c]['length']

		if gc > avg_gcl:
			gc_diff = gc - avg_gcl
		else:
			gc_diff = avg_gcl - gc	
		score = (gc_diff - cov)*length	
		total_score += score		
	return total_score/total_len	

contigs_dict = {}
str_list = read_file(details)
for string in str_list:
	string = string.split("\t")
	if "_" in string[0]:
		name = string[0]
		cov = string[1]
		rd = string[2]
		gc = string[3]
		length = string[4]
		contigs_dict[name] = {}
		contigs_dict[name]['cov'] = float(cov) 
		contigs_dict[name]['rd'] = float(rd) 
		contigs_dict[name]['gc'] = float(gc) 
		contigs_dict[name]['length'] = float(length) 

ILP_dict = {}
str_list = read_file(ILP_plasmids)
for string in str_list:
	string = string.split(";")
	plasmid = string[0]
	ILP_dict[plasmid] = {}
	seq	= string[1].split(",")
	seq = [c[:-1] for c in seq]

	ILP_dict[plasmid]['seq'] = seq  
	ILP_dict[plasmid]['obj'] = eval_obj(seq, contigs_dict)

#print(ILP_dict)

greedy_dict = {}
str_list = read_file(greedy_file)
for string in str_list:
	string = string.split(";")
	plasmid = string[0]
	greedy_dict[plasmid] = {}
	seq	= string[1].split(",")
	seq = [c[:-1]+"_0" for c in seq]

	greedy_dict[plasmid]['seq'] = seq  
	greedy_dict[plasmid]['obj'] = eval_obj(seq, contigs_dict)

real_dict = {}
str_list = read_file(real_file)
for string in str_list:
	string = string.split(";")
	plasmid = string[0]
	real_dict[plasmid] = {}
	seq	= string[1].split(",")
	seq = [c[:-1]+"_0" for c in seq]

	real_dict[plasmid]['seq'] = seq  
	real_dict[plasmid]['obj'] = eval_obj(seq, contigs_dict)	

print("\nILP plasmids")
for plasmid in ILP_dict:
	print(plasmid)
	print(ILP_dict[plasmid]['seq'])
	print(ILP_dict[plasmid]['obj'])
	print("\n")

print("\nGreedy plasmids")
for plasmid in greedy_dict:
	print(plasmid)
	print(greedy_dict[plasmid]['seq'])
	print(greedy_dict[plasmid]['obj'])
	print("\n")

print("\nReal plasmids")
for plasmid in real_dict:
	print(plasmid)
	print(real_dict[plasmid]['seq'])
	print(real_dict[plasmid]['obj'])
	print("\n")


