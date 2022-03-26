__author__ = 'amane'

#-------------------

#USAGE: 
#time python2.7 plasmids_iterative.py sample_id alpha_1 alpha_2 alpha_3 

from gurobipy import *
from sys import argv
import os
import time
from random import randint
import math

import plasmids_preprocessing
import plasmids_postprocessing
import rmcircular

sample_dir = '/home/aniket/PhD/Plasmids-Optimization/data/unicycler_pipeline/'
output_dir = '/home/aniket/PhD/Plasmids-Optimization/exp/2021-06-14__iterative_MILP_edge_novelty/output/' 

#sample_dir = '/home/aniket/PhD/Plasmids-Optimization/exp/2021-01-29__iterative_MILP_modified_graph/test_examples/' 
#output_dir = '/home/aniket/PhD/Plasmids-Optimization/exp/2021-01-29__iterative_MILP_modified_graph/test_output/' 


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
output_folder = output_dir + 'sample_' + sample_id + '/' + ratios + '/nplasmids_'+str(nplasmids) + '/MILP'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
#output_filename = argv[6]
output_filename = 'putative_plasmids.fasta'
questionable_filename = 'questionable_plasmids.fasta'
score_filename = 'MILP_objective.csv'
contigs_filename = 'contig_chains.csv'
q_contigs_filename = 'questionable_contig_chains.csv'

links_filename = 'links.csv'

#-----------------------------------------------
#Main program
contigs_dict = {}
links_list = []
seeds_set = set()
rd_dict = {}
edge_fam_dict = {}
contigs_dict, links_list, rd_dict, edge_fam_dict = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list, rd_dict, edge_fam_dict)
seeds_set = plasmids_preprocessing.get_seeds(seeds_file, seeds_set, rd_dict)




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

contigs_dict = plasmids_preprocessing.get_gene_coverage(mapping_file, contigs_dict, rd_dict)
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
	#print(c, contigs_dict[c]['Gene_coverage'], contigs_dict[c]['Read_depth'], contigs_dict[c]['Length'])	
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
links_file = open(os.path.join(output_folder, links_filename), "w")

logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"w")

n_iter = 0
q_iter = 0


#print(extr_dict)
#print(links_list)
iter_count = 0
while len(seeds_set) > 0:
#while iter_count <= 1:
	iter_count += 1
	print("\n\n\n\n\n")
	#-----------------------------------------------
	#Initializing the ILP
	m = Model("Plasmids")
	m.params.LogFile= os.path.join(output_folder,'m.log')
	m.setParam(GRB.Param.TimeLimit, 240.0)
	#m.setParam('MIPGap', 0.05)
	m.setParam(GRB.Param.MIPGap, 0.05)
	#Minimum seed read depth check
	#-----------------------------------------------
	#for x in seeds_set:
	#	print(x, contigs_dict[x]['Read_depth'])




	contigs = {}
	degree = {}
	contigs_ext = {}
	ceil = {}
	contigs, degree, contigs_ext, ceil = plasmids_preprocessing.contig_vars(m, contigs_dict, contigs, degree, contigs_ext, ceil, nplasmids)

	links = {}
	links = plasmids_preprocessing.link_vars(m, links_list, links, nplasmids)

	rd_diff, counted_rd_diff = {}, {}	
	rd_diff, counted_rd_diff \
		= plasmids_preprocessing.rd_vars(m, contigs_dict, rd_diff, counted_rd_diff, nplasmids)

	mean_GC, counted_GC_mean, GC_diff, counted_GC_diff = {}, {}, {}, {}
	mean_GC, counted_GC_mean, GC_diff, counted_GC_diff \
		= plasmids_preprocessing.GC_vars(m, contigs_dict, mean_GC, counted_GC_mean, GC_diff, counted_GC_diff, nplasmids)

	counted_ln = {}
	counted_ln = plasmids_preprocessing.ln_vars(m, contigs_dict, counted_ln, nplasmids)

	counted_seed = {}
	counted_seed = plasmids_preprocessing.seed_vars(m, contigs_dict, counted_seed, nplasmids)

	#wtd_gd, counted_wtd_gd = {}, {}
	#wtd_gd, counted_wtd_gd \
	#	= plasmids_preprocessing.gd_vars(m, contigs_dict, wtd_gd, counted_wtd_gd, nplasmids)

	#print("\nLinks list:")
	#for link in links[0]:
	#	print(link)	
	#print("\n")	

	#-----------------------------------------------
	#Setting up the expression for the objective function
	expr = LinExpr()
	alpha1, alpha2, alpha3 = float(alpha1), float(alpha2), float(alpha3) 

	for p in rd_diff:
		for c in rd_diff[p]:
			#expr.addTerms(alpha1, counted_rd_diff[p][c])

			#expr.addTerms(alpha2, contigs[p][c])

			expr.addTerms(-alpha2*contigs_dict[c]['Gene_coverage'], contigs[p][c])

			expr.addTerms(alpha3, counted_GC_diff[p][c])
		#expr.addConstant(2)	


	m.setObjective(expr, GRB.MINIMIZE)

	#-----------------------------------------------
	#Setting up constraints

	constraint_count = 0

	#Weeding out zero gene density contigs at the end
	#for p in degree:
		#print(contigs_ext[p])
	#	for c in degree[p]:
	#		end1 = (c, 'h')
	#		end2 = (c, 't')
	#		print(c, end1, end2)			
	#		m.addConstr(degree[p][c] >= contigs_ext[p][end1] - contigs_ext[p][ext2], "degree_lbd1")
	#		m.addConstr(degree[p][c] >= contigs_ext[p][end2] - contigs_ext[p][ext1], "degree_lbd2")
	#		m.addConstr(degree[p][c] + 1 - contigs_dict[c]['Density'] <= 1.9, "end_contig_gd")

	#Constraint type 1
	#A link 'e' is in 'p' only if both its endpoints are in 'p'
	for p in links:
		for e in links[p]:
			end1, end2 = e[0], e[1]
			m.addConstr(links[p][e] <= contigs_ext[p][end1], "link_ubd")
			m.addConstr(links[p][e] <= contigs_ext[p][end2], "link_ubd")
			constraint_count += 2

			

	#Constraint type 2
	#Consistency: A contig extremity can occur only once in a plasmid
	for p in contigs_ext:
		for extr in contigs_ext[p]:
			c = extr[0]
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
			link_count = 0
			for link in extr_dict[extr]:
				if link in links[p]:
					link_count += 1
					expr.addTerms(1, links[p][link])
			#if link_count > 0:
			#	expr = expr/link_count
			m.addConstr(contigs_ext[p][extr] <= expr + 0.99, "extr_ubd")
			m.addConstr(contigs_ext[p][extr] >= expr, "extr_lbd")
			constraint_count += 2

	#Constraint type 4
	#A contig 'c' is in 'p' if at least one of its endpoints is in 'p'
	for p in contigs:
		for c in contigs[p]:
			end1, end2 = (c, 'h'), (c, 't')
			m.addConstr(contigs[p][c] >= contigs_ext[p][end1], "contig_lbd1")
			m.addConstr(contigs[p][c] >= contigs_ext[p][end2], "contig_lbd2")
			constraint_count += 2						
	
	
	#Constraint type 5
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
		m.addConstr(expr <= 175000, "plasmid_length_upper_bound")
		constraint_count += 1	

	#Constraint type 6
	#counted-seed = seed * existence-of-a-contig-in-plasmid
	for p in counted_seed:
		for c in counted_seed[p]:
			if c in contigs[p]:
				m.addConstr(counted_seed[p][c] == contigs_dict[c]['Seed']*contigs[p][c], "counted_seed")
				constraint_count += 1				

	#Constraint type 8
	#Diff value of average read depth and plasmid read-depth
	for p in rd_diff:
		for c in rd_diff[p]:
			if c in rd_diff[p]:
				m.addConstr(rd_diff[p][c] >= contigs[p][c] - contigs_dict[c]['Read_depth']/rd_graph, "rd-diff_lbd1")
				m.addConstr(rd_diff[p][c] >= 0, "rd-diff_lbd2")
				constraint_count += 2

	for p in counted_rd_diff:
		for c in counted_rd_diff[p]:
			expr = LinExpr()
			expr.addTerms(0.01, contigs[p][c])
			m.addConstr(ceil[p][c] <= expr + 0.99, "get-ceil-ubd")
			m.addConstr(ceil[p][c] >= expr, "get-ceil-lbd")

			m.addConstr(counted_rd_diff[p][c] <= UBD_rd*ceil[p][c], "counted-rd-diff1")
			m.addConstr(counted_rd_diff[p][c] <= rd_diff[p][c], "counted-rd-diff2")
			m.addConstr(counted_rd_diff[p][c] >= rd_diff[p][c] - (1 - ceil[p][c])*UBD_rd, "counted-rd-diff3")
			m.addConstr(counted_rd_diff[p][c] >= 0, "counted-rd-diff4")
			constraint_count += 4

	#Constraint type 9
	#Diff value of computed mean GC cont and plasmid GC cont
	for p in counted_GC_mean:
		expr1 = LinExpr()
		expr2 = LinExpr()		
		for c in counted_GC_mean[p]:
			expr1.addTerms(contigs_dict[c]['Length'], counted_GC_mean[p][c])
			expr2.addTerms(contigs_dict[c]['Length']*contigs_dict[c]['GC_cont'], contigs[p][c])
		m.addConstr(expr1 == expr2, "counted-GC-mean")
		constraint_count += 1			

	for p in counted_GC_mean:
		for c in counted_GC_mean[p]:
			m.addConstr(counted_GC_mean[p][c] <= UBD_GC*contigs[p][c], "ctd-GC-mean1")
			m.addConstr(counted_GC_mean[p][c] <= mean_GC[p], "ctd-GC-mean2")
			m.addConstr(counted_GC_mean[p][c] >= mean_GC[p] - UBD_GC*(1-contigs[p][c]), "ctd-GC-mean3")
			m.addConstr(counted_GC_mean[p][c] >= 0, "ctd-GC-mean4")	
			constraint_count += 4	

	for p in GC_diff:
		for c in GC_diff[p]:
			m.addConstr(GC_diff[p][c] >= mean_GC[p] - contigs_dict[c]['GC_cont']*contigs[p][c], "GC-diff_lbd1")
			m.addConstr(GC_diff[p][c] >= contigs_dict[c]['GC_cont']*contigs[p][c] - mean_GC[p], "GC-diff_lbd2")
			constraint_count += 2

	for p in counted_GC_diff:
		for c in counted_GC_diff[p]:
			m.addConstr(counted_GC_diff[p][c] <= UBD_GC*contigs[p][c], "counted-GC-diff1")
			m.addConstr(counted_GC_diff[p][c] <= GC_diff[p][c], "counted-GC-diff2")
			m.addConstr(counted_GC_diff[p][c] >= GC_diff[p][c] - (1 - contigs[p][c])*UBD_GC, "counted-GC-diff3")
			m.addConstr(counted_GC_diff[p][c] >= 0, "counted-GC-diff4")
			constraint_count += 4					
				

	#Constraint type 10
	#Each plasmids should have at least one seed
	for p in counted_seed:
		expr = LinExpr()
		for c in counted_seed[p]:
			expr.addTerms(1, counted_seed[p][c])
		m.addConstr(expr >= 1, "seed_existence")
		constraint_count += 1


	#Constraint type 11
	#Preliminary path constraints
	for p in contigs:
		c_expr = LinExpr()
		for c in contigs[p]:
			c_expr.addTerms(1, contigs[p][c])
		l_expr = LinExpr()
		for e in links[p]:
			l_expr.addTerms(1, links[p][e])
		m.addConstr(c_expr == l_expr + 1, "path_simple")

	#Constraint type 12
	#Edge novelty constraints
	for p in links:
		for edge in edge_fam_dict:
			edge_expr = LinExpr()
			for e in edge_fam_dict[edge]:
				if e in links[p]:
					edge_expr.addTerms(1, links[p][e])
			m.addConstr(edge_expr <= 1, "edge_novelty")
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
		
		m.addConstr(mean_GC[p] >= 0, "non-negativity")

		for c in counted_ln[p]:
			m.addConstr(counted_ln[p][c] >= 0, "non-negativity")
			m.addConstr(counted_seed[p][c] >= 0, "non-negativity")
			m.addConstr(rd_diff[p][c] >= 0, "non-negativity")
			m.addConstr(GC_diff[p][c] >= 0, "non-negativity")





	#while len(seeds_set) > 0:
	#if m.status == GRB.Status.INF_OR_UNBD:
		# Turn presolve off to determine whether model is infeasible
		# or unbounded
	#	print("It is not optimal")
	#	m.setParam(GRB.Param.Presolve, 0)
	#	m.optimize()

	#if m.status == GRB.Status.OPTIMAL:
	#status = m.status
	logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"a")
	circular = 1
	i = 1
	while circular == 1:
		start = time.time()
		m.optimize()
		stop = time.time()
		duration = stop - start
		#time_file.write("Iteration "+str(i)+":\t"+str(duration)+"\n")

		if m.status == GRB.Status.INFEASIBLE:
			print ('The model cannot be solved because it is infeasible')
		elif m.status == GRB.Status.UNBOUNDED:
			print ('The model cannot be solved because it is unbounded')
		elif m.status == GRB.Status.INF_OR_UNBD:
			print ('The model cannot be solved because it is infeasible or unbounded ')

		if m.status == GRB.Status.INF_OR_UNBD or m.status == GRB.Status.INFEASIBLE or m.status == GRB.Status.UNBOUNDED:
			print ('The model cannot be solved because it is infeasible or unbounded ')
			m.computeIIS()
			m.write("m.ilp")
			for con in m.getConstrs():
				if con.IISConstr:
					print('%s' % con.constrName) 
			exit (1)

		for p in contigs:
			for c in contigs[p]:
				if contigs[p][c].x > 0:	
					print("Before removing cycles: ", c, contigs[p][c].x)
			for e in links[p]:
				if links[p][e].x > 0:
					print("Before removing cycles: ", e, links[p][e].x)			

		print("\n")	

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
					if end1 not in soln_ext_dict[p]:
						soln_ext_dict[p][end1] = [end2]
					else:
						soln_ext_dict[p][end1].append(end2)
					if end2 not in soln_ext_dict[p]:
						soln_ext_dict[p][end2] = [end1]
					else:
						soln_ext_dict[p][end2].append(end1)	

			solution_seq[p], contig_chain, plasmid_parts = plasmids_postprocessing.get_seq(solution_links[p], soln_ext_dict[p], contigs_dict, contigs)			

		circ_seq = {}
		k = 1 
		logfile_3.write("Iteration "+ str(i)+"\n")

		#logfile_4 = open(os.path.join(output_folder,'Final_soln.log'),"w")
		circular = 0	

		circular, circ_seq, solution_seq[p] = \
			rmcircular.check_if_circular(solution_seq, circ_seq, k, i, soln_ext_dict, solution_links, contigs, logfile_3)
		
		if circular == 1:	#If chr is circular, adds corresponding constraint
			#x = randint(1, len(circ_seq))

			for x in circ_seq:

				chosen_seq = set()
				#for e in circ_seq[x]['Seq']:
				#	if e not in chosen_seq:
				#		chosen_seq.add(e)
				#p = circ_seq[x]['Plasmid']		
				chosen_seq = circ_seq[x]['Seq']
				print("Chosen circ sequence: ", chosen_seq)
				logfile_3.write("\nChosen seed at iteration "+str(i)+": "+str(x))
				logfile_3.write("\nChosen chr at iteration "+str(i)+"\n")
				expr = LinExpr()
				for e in chosen_seq:
	
					print("Simple",e)
					if e in links[p]:
						expr.addTerms(1, links[p][e])
					elif e[::-1] in links[p]:
						expr.addTerms(1, links[p][e[::-1]])
				logfile_3.write("\n\n")		
				m.addConstr(expr <= len(chosen_seq) - 1, "circular")
						
				print("\n", i, "th circular constraint added\n")
				i += 1
			#else:
			#	circular = 0
				

		if i >= 5:
			break	
#_______________

	print('Obj:', m.objVal)

	plasmid_length = {}
	plasmid_gd = {}
	plasmid_rd = {}
	for p in contigs:
		plasmid_length[p] = 0
		plasmid_gd[p] = 0
		plasmid_rd[p] = 0
		for c in contigs[p]:
			if contigs[p][c].x > 0:
				print(contigs[p][c].varName, contigs[p][c].x)
				#print(counted_rd[p][c].varName, counted_rd[p][c].x)
				print(contigs_dict[c]['Read_depth'], contigs_dict[c]['Length'], contigs_dict[c]['Gene_coverage'])
				plasmid_length[p] += contigs_dict[c]['Length']
				plasmid_gd[p] += contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length']
				if contigs_dict[c]['Read_depth'] <= rd_graph:
					plasmid_rd[p] += contigs_dict[c]['Read_depth']*contigs_dict[c]['Length']
				else:
					plasmid_rd[p] += rd_graph*contigs_dict[c]['Length']
		plasmid_gd[p] = plasmid_gd[p]/plasmid_length[p]
		plasmid_rd[p] = plasmid_rd[p]/plasmid_length[p]
		for e in links[p]:
			if links[p][e].x > 0:
				print(links[p][e].varName, links[p][e].x)		
		#print("Plasmid length: ", p, plasmid_length[p])
		#print("Plasmid gene density: ", p, plasmid_gd[p])
		#print("Plasmid depth of coverage: ", p, plasmid_rd[p])
		#print("\n")



	output_file = open(os.path.join(output_folder, output_filename), "a")
	questionable_file = open(os.path.join(output_folder, questionable_filename), "a")
	score_file = open(os.path.join(output_folder, score_filename), "a")
	contigs_file = open(os.path.join(output_folder, contigs_filename), "a")
	q_contigs_file = open(os.path.join(output_folder, q_contigs_filename), "a")

	links_file = open(os.path.join(output_folder, links_filename), "a")

	solution_links = {}		
	solution_seq = {}
	soln_ext_dict = {}
	for p in links:
		solution_links[p] = set()
		solution_seq[p] = ''
		soln_ext_dict[p] = {}
		for e in links[p]:
			if links[p][e].x > 0.5:
		
				end1, end2 = e[0], e[1]
				c1, ext1 = end1[0], end1[1]
				c2, ext2 = end2[0], end2[1]			

				solution_links[p].add(e)
				if end1 not in soln_ext_dict[p]:
					soln_ext_dict[p][end1] = [end2]
				else:
					soln_ext_dict[p][end1].append(end2)
				if end2 not in soln_ext_dict[p]:
					soln_ext_dict[p][end2] = [end1]
				else:
					soln_ext_dict[p][end2].append(end1)	

		#Getting solution seq and contig chain			
		solution_seq[p], contig_chain, plasmid_parts = plasmids_postprocessing.get_seq(solution_links[p], soln_ext_dict[p], contigs_dict, contigs)
		#Contig chain for single-contig plasmids
		if len(contig_chain) == 0:
			for c in contigs[p]:
				if contigs[p][c].x > 0:
					contig_chain.append(str(c)+'+')
		print("Selected chain: ", contig_chain)

		contig_chain, plasmid_length[p], plasmid_gd[p], plasmid_rd[p], solution_seq[p] = plasmids_postprocessing.refine_chain(contig_chain, contigs_dict)	

		#Sorting into putative and questionable
		if plasmid_length[p] >= 1500 and plasmid_gd[p] >= 0.3:
			output_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length[p])+"\t"+"gene_density="+str(plasmid_gd[p])+"\t"+"mean_read_depth="+str(plasmid_rd[p])+"\n")
			if len(solution_seq[p]) != 0:							#if length of contig chain is more than zero
				output_file.write(solution_seq[p]+"\n")
			else:													#if plasmid is a single contig
				c = contig_chain[0][:-1]								
				output_file.write(contigs_dict[c]['Sequence']+"\n")
			
			#Printing putative contig chains to file
			contigs_file.write("plasmid_"+str(n_iter)+";")
			links_file.write("plasmid_"+str(n_iter)+";\n")
			for x in contig_chain:
				if x == contig_chain[0]:
					contigs_file.write(x)
				else:
					contigs_file.write(","+x)					
			contigs_file.write("\n")
			if len(contig_chain) == 1:
				links_file.write(contig_chain[0])
				links_file.write("\n")
			#Printing putative contig chain details to file	
			contigs_file.write("#Contig\tGene_coverage\tGC_cont\tRead_depth\tLength\n")		
			for x in contig_chain:
				c = x[:-1]
				contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['Gene_coverage'])+"\t"+str(contigs_dict[c]['GC_cont'])+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
			contigs_file.write("\n")

			
			for x in plasmid_parts:
				for i in range(len(x)):
					y = x[i]
					if i == 0:
						links_file.write(str(y[0]+y[1]))
					else:
						links_file.write(","+str(y[0]+y[1]))	
				links_file.write("\n")			
				print(x, "PLASMID PART")

			#Recording scores
			rd_sum, gd_sum, GC_sum = 0, 0, 0
			for p in contigs:
				for c in contigs[p]:
					rd_sum += alpha1*counted_rd_diff[p][c].x
					#expr.addTerms(alpha1*contigs_dict[c]['Length'], counted_diff[p][c])
					gd_sum += alpha2*(-contigs_dict[c]['Gene_coverage'])*contigs[p][c].x
					#print(c, counted_wtd_gd[p][c].x)
					#if contigs[p][c].x > 0:
					#	print("Selected: ",c, counted_wtd_gd[p][c].x, counted_wtd_GC_diff[p][c].x, counted_wtd_rd_diff[p][c].x)
					#expr.addTerms(-alpha2*contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length'], counted_rd[p][c])
					#GC_sum += alpha3*counted_GC_diff[p][c].x
					#expr.addTerms(-alpha3*(GC_mean-contigs_dict[c]['GC_cont']), counted_ln[p][c])
			score_file.write("putative_plasmid_"+str(n_iter)+"\t\t"+str(m.objVal)+"\t"+str(rd_sum)+"\t"+str(gd_sum)+"\t"+str(GC_sum)+"\n")
			#score_file.write("plasmid_"+str(n_iter)+"\t"+str(m.objVal)+"\t"+str(wtd_rd_diff[p].x)+"\t"+str(wtd_gd[p].x)+"\t"+str(wtd_GC_diff[p].x)+"\n")
			n_iter += 1


		
		else:
			#Same for questionable plasmids		
			questionable_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length[p])+"\t"+"gene_density="+str(plasmid_gd[p])+"\t"+"mean_read_depth="+str(plasmid_rd[p])+"\n")
			if len(solution_seq[p]) != 0:
				questionable_file.write(solution_seq[p]+"\n")
			else:
				c = contig_chain[0][:-1]
				questionable_file.write(contigs_dict[c]['Sequence']+"\n")
	
			q_contigs_file.write("plasmid_"+str(q_iter)+";")
			for x in contig_chain:
				if x == contig_chain[0]:
					q_contigs_file.write(x)

				else:
					q_contigs_file.write(","+x)
			q_contigs_file.write("\n")		
			for x in contig_chain:
				c = x[:-1]
				q_contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['Gene_coverage'])+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
			q_contigs_file.write("\n")

			rd_sum, gd_sum, GC_sum = 0, 0, 0
			for p in contigs:
				for c in contigs[p]:
					rd_sum += alpha1*counted_rd_diff[p][c].x
					#expr.addTerms(alpha1*contigs_dict[c]['Length'], counted_diff[p][c])
					gd_sum += alpha2*(-contigs_dict[c]['Gene_coverage'])*contigs[p][c].x
					#print(c, counted_wtd_gd[p][c].x)
					#if contigs[p][c].x > 0:
					#	print("Selected: ",c, counted_wtd_gd[p][c].x, counted_wtd_GC_diff[p][c].x, counted_wtd_rd_diff[p][c].x)
					#expr.addTerms(-alpha2*contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length'], counted_rd[p][c])
					#GC_sum += alpha3*counted_GC_diff[p][c].x
					#expr.addTerms(-alpha3*(GC_mean-contigs_dict[c]['GC_cont']), counted_ln[p][c])
			score_file.write("questionable_plasmid_"+str(n_iter)+"\t"+str(m.objVal)+"\t"+str(rd_sum)+"\t"+str(gd_sum)+"\t"+str(GC_sum)+"\n")
			#score_file.write("plasmid_"+str(n_iter)+"\t"+str(m.objVal)+"\t"+str(wtd_rd_diff[p].x)+"\t"+str(wtd_gd[p].x)+"\t"+str(wtd_GC_diff[p].x)+"\n")
			q_iter += 1
	

	#Updating assembly graph and formulation
	for p in contigs:
		#print("Mean GC cont: ", mean_GC[p].x)
		for c in contigs[p]:
			if contigs[p][c].x > 0:
				#print(c, counted_GC_mean[p][c].x)
				#print(c, GC_diff[p][c].x)
				#print(c, counted_GC_diff[p][c].x)

				print(c,rd_diff[p][c].x)
				print(c,counted_rd_diff[p][c].x)
				print(c, ceil[p][c].x)

				#print(c, 'h', contigs_ext[p][(c, 'h')].x)
				#print(c, 't', contigs_ext[p][(c, 't')].x)

				print("BEFORE UPDATING: ", c, contigs_dict[c]['Read_depth'])
				contigs_dict[c]['Read_depth'] = max(0, contigs_dict[c]['Read_depth'] - rd_graph*contigs[p][c].x)
				print("AFTER UPDATING: ", c, contigs_dict[c]['Read_depth'])

				if c in seeds_set and contigs_dict[c]['Read_depth'] <= 0.5:
					contigs_dict[c]['Seed'] = 0
					seeds_set.remove(c)

				if contigs_dict[c]['Read_depth'] <= 0.05:
					print("Deleting contig ", c)
					del(contigs_dict[c])				

					temp_list = []
					for e in links_list:						
						if (c, 'h') not in e and (c, 't') not in e:
							#print(e)
							temp_list.append(e)
					links_list = temp_list


						


