__author__ = 'amane'

#-------------------
#The program takes as input, an assembly graph and a file mapping genes to contigs and outputs a set of plasmids in the form of sets of contigs.
#The number of plasmids "nplasmids" is to be provided as part of the input.

#USAGE: 
#time python2.7 plasmids_main.py nplasmids samples/sample_name/assembly.gfa 
#samples/sample_name/filtered_genes_to_contigs.csv samples/sample_name/seed_contigs.csv 
#Here, nplasmids is the number of expected plasmids

from gurobipy import *
from sys import argv
import time

import plasmids_preprocessing
import plasmids_postprocessing

#-----------------------------------------------
#Number of plasmids
nplasmids = int(argv[1]) 

#-----------------------------------------------
#Input files (Assembly graph, gene-contig mapping, seed contigs)
assembly_file = argv[2]
mapping_file = argv[3]
seeds_file = argv[4]

#-----------------------------------------------
#Output file 
output_folder = argv[5]
output_folder = output_folder+'_'+str(nplasmids)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
output_filename = argv[6]

logfile = "MILP.log"

#-----------------------------------------------
#Main program
contigs_dict = {}
links_list = []
seeds_set = set()

contigs_dict, links_list = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list)
seeds_set = plasmids_preprocessing.get_seeds(seeds_file, seeds_set)

GC_total = 0
ln_total = 0
for c in contigs_dict:
	if c in seeds_set:
		contigs_dict[c]['Seed'] = 1
	else:
		contigs_dict[c]['Seed'] = 0

	GC_total += contigs_dict[c]['GC_cont']*contigs_dict[c]['Length']
	ln_total += contigs_dict[c]['Length']
	
GC_mean = GC_total/ln_total		

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

contigs_dict = plasmids_preprocessing.get_gene_coverage(mapping_file, contigs_dict)

#-----------------------------------------------
#Initializing the ILP
m = Model("Plasmids")

contigs = {}
contigs_ext = {}
contigs, contigs_ext = plasmids_preprocessing.contig_vars(m, contigs_dict, contigs, contigs_ext, nplasmids)

links = {}
links = plasmids_preprocessing.link_vars(m, links_list, links, nplasmids)

rd = {}
counted_rd = {}
mean_rd = {}
counted_mean_rd = {}
rd, counted_rd, mean_rd, counted_mean_rd \
	= plasmids_preprocessing.rd_vars(m, contigs_dict, rd, counted_rd, mean_rd, counted_mean_rd, nplasmids)

counted_ln = {}
counted_ln = plasmids_preprocessing.ln_vars(m, contigs_dict, counted_ln, nplasmids)

counted_seed = {}
counted_seed = plasmids_preprocessing.seed_vars(m, contigs_dict, counted_seed, nplasmids)

diff, counted_diff = {}, {}
diff, counted_diff = plasmids_preprocessing.diff_vars(m, contigs_dict, diff, counted_diff, nplasmids)


#-----------------------------------------------
#Setting up the expression for the objective function
expr = LinExpr()
#alpha1, alpha2, alpha3 = 0.25, 0.25, 0.25
#alpha1, alpha2, alpha3 = 0.25, 0, 0
#alpha1, alpha2, alpha3 = 0, 0.25, 0
#alpha1, alpha2, alpha3 = 0, 0, 0.25
#alpha1, alpha2, alpha3 = 0.25, 0.25, 0
#alpha1, alpha2, alpha3 = 0.25, 0, 0.25
#alpha1, alpha2, alpha3 = 0, 0.25, 0.25
#alpha1, alpha2, alpha3 = 0.5, 0.25, 0.25
#alpha1, alpha2, alpha3 = 0.25, 0.5, 0.25
#alpha1, alpha2, alpha3 = 0.25, 0.25, 0.5
#alpha1, alpha2, alpha3 = 2.5, 1.25, 0.25
#alpha1, alpha2, alpha3 = 1.25, 0.5, 0.25
alpha1, alpha2, alpha3 = 0.5, 0.25, 0




for p in diff:
	for c in diff[p]:
		expr.addTerms(alpha1*contigs_dict[c]['Length'], counted_diff[p][c])

		expr.addTerms(-alpha2*contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length'], counted_rd[p][c])

		expr.addTerms(-alpha3*(GC_mean-contigs_dict[c]['GC_cont']), counted_ln[p][c])

m.setObjective(expr, GRB.MINIMIZE)

#-----------------------------------------------
#Setting up constraints

constraint_count = 0

#Constraint type 1
#A link 'e' is in 'p' only if both its endpoints are in 'p'
for p in links:
	for e in links[p]:
		end1, end2 = e[0], e[1]
		if end1 in contigs_ext[p] and end2 in contigs_ext[p]:
			m.addConstr(links[p][e] >= (contigs_ext[p][end1] + contigs_ext[p][end2])/2 - 0.99, "link_lbd")
			m.addConstr(links[p][e] <= (contigs_ext[p][end1] + contigs_ext[p][end2])/2, "link_ubd")
			constraint_count += 2

#Constraint type 2
#Consistency: A contig extremity can occur only once in a plasmid
for p in contigs_ext:
	for extr in contigs_ext[p]:
		expr = LinExpr()
		for link in extr_dict[extr]:
			if link in links[p]:
				expr.addTerms(1, links[p][link])
		m.addConstr(expr <= 1, "consistency")
		constraint_count += 1

#Constraint type 3
#An extremity is in 'p' only if at least one edge is incident on it
for p in contigs_ext:
	for extr in contigs_ext[p]:
		expr = LinExpr()
		for link in extr_dict[extr]:
			if link in links[p]:
				expr.addTerms(1, links[p][link])
		if len(extr_dict[extr]) > 0:
			expr = expr/len(extr_dict[extr])
		m.addConstr(contigs_ext[p][extr] <= expr + 0.99, "extr_ubd")
		m.addConstr(contigs_ext[p][extr] >= expr, "extr_lbd")
		constraint_count += 2

#Constraint type 4
#A contig 'c' is in 'p' if at least one of its endpoints is in 'p'
for p in contigs:
	for c in contigs[p]:
		end1, end2 = (c, 'h'), (c, 't')
		m.addConstr(contigs[p][c] <= (contigs_ext[p][end1] + contigs_ext[p][end2])/2 + 0.99, "contig_ubd")
		m.addConstr(contigs[p][c] >= (contigs_ext[p][end1] + contigs_ext[p][end2])/2, "contig_lbd")
		constraint_count += 2			


#Constraint type 5
#Total of read depth contribution to each plasmid is less than contig read depth
for c in contigs_dict:
	expr = LinExpr()
	for p in rd:
		expr.addTerms(1, rd[p][c])
	m.addConstr(expr <= contigs_dict[c]['Read_depth'], "read-depth_cap")
	constraint_count += 1		

#Constraint type 6
#counted-read-depth = read-depth * existence-of-a-contig-in-plasmid
for p in counted_rd:
	for c in counted_rd[p]:
		if c in contigs[p]:
			m.addConstr(counted_rd[p][c] <= contigs_dict[c]['Read_depth']*contigs[p][c], "counted_read-depth_1")
			m.addConstr(counted_rd[p][c] <= rd[p][c], "counted_read-depth_2")
			m.addConstr(counted_rd[p][c] >= rd[p][c] - (1 - contigs[p][c])*contigs_dict[c]['Read_depth'], "counted_read-depth_3")
			m.addConstr(counted_rd[p][c] >= 0, "counted_read-depth_4")
			constraint_count += 4				

#Constraint type 7
#counted-length = length * existence-of-a-contig-in-plasmid
for p in counted_ln:
	for c in counted_ln[p]:
		if c in contigs[p]:
			m.addConstr(counted_ln[p][c] == contigs_dict[c]['Length']*contigs[p][c], "counted_length")
			constraint_count += 1

for p in contigs:
	expr = LinExpr()
	for c in contigs[p]:
		expr.addTerms(1, counted_ln[p][c])
	m.addConstr(expr <= 200000, "plasmid_length_upper_bound")	

#Constraint type 8
#counted-seed = seed * existence-of-a-contig-in-plasmid
for p in counted_seed:
	for c in counted_seed[p]:
		if c in contigs[p]:
			m.addConstr(counted_seed[p][c] == contigs_dict[c]['Seed']*contigs[p][c], "counted_seed")
			constraint_count += 1			

#Constraint type 9
#Computing mean read depth for each plasmid
for p in counted_mean_rd:
	expr1 = LinExpr()
	expr2 = LinExpr()
	for c in counted_mean_rd[p]:
		if c in contigs[p]:
			m.addConstr(counted_mean_rd[p][c] <= 1000000*contigs[p][c], "counted_mean-rd1")
			m.addConstr(counted_mean_rd[p][c] <= mean_rd[p], "counted_mean-rd2")
			m.addConstr(counted_mean_rd[p][c] >= mean_rd[p] - (1 - contigs[p][c])*1000000, "counted_mean-rd3")
			m.addConstr(counted_mean_rd[p][c] >= 0, "counted_mean-rd4")

			expr1.addTerms(contigs_dict[c]['Length'], counted_mean_rd[p][c])
			expr2.addTerms(contigs_dict[c]['Length'], counted_rd[p][c])
			m.addConstr(expr1 == expr2, "computing-mean-rd")
			constraint_count += 5	


#Constraint type 10
#Absolute value (diff) of mean read depth and plasmid read-depth
for p in counted_rd:
	for c in counted_rd[p]:
		if c in diff[p]:
			m.addConstr(diff[p][c] >= mean_rd[p] - counted_rd[p][c], "diff_lbd1")
			m.addConstr(diff[p][c] >= counted_rd[p][c] - mean_rd[p], "diff_lbd2")
			constraint_count += 2

for p in counted_diff:
	for c in counted_diff[p]:
		m.addConstr(counted_diff[p][c] <= 10*contigs[p][c], "counted_diff1")
		m.addConstr(counted_diff[p][c] <= diff[p][c], "counted_diff2")
		m.addConstr(counted_diff[p][c] >= diff[p][c] - (1 - contigs[p][c])*10, "counted_diff3")
		m.addConstr(counted_diff[p][c] >= 0, "counted_diff4")


#Constraint type 11
#Each plasmids should have at least one seed
for p in counted_seed:
	expr = LinExpr()
	for c in counted_seed[p]:
		expr.addTerms(1, counted_seed[p][c])
	m.addConstr(expr >= 1, "seed_existence")
	constraint_count += 1


#Constraint type 0
#Non-negativity constraints
for p in contigs:
	for c in contigs[p]:
		m.addConstr(contigs[p][c] >= 0, "non-negativity")
	for extr in contigs_ext[p]:
		m.addConstr(contigs_ext[p][extr] >= 0, "non-negativity")
	for e in links[p]:
		m.addConstr(links[p][e] >= 0, "non-negativity")
	for c in rd[p]:
		m.addConstr(rd[p][c] >= 0, "non-negativity")
		m.addConstr(counted_rd[p][c] >= 0, "non-negativity")
		m.addConstr(counted_mean_rd[p][c] >= 0, "non-negativity")
	m.addConstr(mean_rd[p] >= 0, "non-negativity")
	for c in counted_ln[p]:
		m.addConstr(counted_ln[p][c] >= 0, "non-negativity")
		m.addConstr(counted_seed[p][c] >= 0, "non-negativity")
		m.addConstr(diff[p][c] >= 0, "non-negativity")


m.optimize()

#m.computeIIS()
#m.write("m.ilp")	

plasmid_length = {}
for p in contigs:
	plasmid_length[p] = 0
	for c in contigs[p]:
		if contigs[p][c].x > 0:
			print(contigs[p][c].varName, contigs[p][c].x)
			print(counted_rd[p][c].varName, counted_rd[p][c].x)
			plasmid_length[p] += contigs_dict[c]['Length']
	print(p, plasmid_length[p])
	print(mean_rd[p].varName, mean_rd[p].x)
	print("\n")	


#output_filename = output_filename + "_" + str(nplasmids) + ".fasta"
output_file = open(os.path.join(output_folder, output_filename), "w")

solution_links = {}		
solution_seq = {}
soln_ext_dict = {}
for p in links:
	solution_links[p] = set()
	solution_seq[p] = ''
	soln_ext_dict[p] = {}
	for e in links[p]:
		if links[p][e].x > 0:
			end1, end2 = e[0], e[1]
			c1, ext1 = end1[0], end1[1]
			c2, ext2 = end2[0], end2[1]			

			solution_links[p].add(e)
			soln_ext_dict[p][e[0]] = e[1]
			soln_ext_dict[p][e[1]] = e[0]
	solution_seq[p] = plasmids_postprocessing.get_seq(solution_links[p], soln_ext_dict[p], contigs_dict)
	output_file.write(">plasmid_"+str(p)+"\t"+"length="+str(plasmid_length[p])+"\t"+"mean_read_depth="+str(mean_rd[p].x)+"\n")
	output_file.write(solution_seq[p]+"\n")

	