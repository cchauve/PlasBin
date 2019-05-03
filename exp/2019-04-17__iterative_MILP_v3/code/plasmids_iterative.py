__author__ = 'amane'

#-------------------

#USAGE: 
#time python2.7 plasmids_main.py sample_id alpha_1 alpha_2 alpha_3 

from gurobipy import *
from sys import argv
import os
import time
from random import randint

import plasmids_preprocessing
import plasmids_postprocessing
import rmcircular

sample_dir = '/home/aniket/python_scripts/Plasmids/data/unicycler_pipeline/'
#sample_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2019-04-17__iterative_MILP_v3/'
output_dir = '/home/aniket/python_scripts/Plasmids-Optimization/exp/2019-04-17__iterative_MILP_v3/output/' 

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

#-----------------------------------------------
#Main program
contigs_dict = {}
links_list = []
seeds_set = set()

contigs_dict, links_list = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list)
seeds_set = plasmids_preprocessing.get_seeds(seeds_file, seeds_set)

#print(seeds_set)

#GC_total = 0
ln_total = 0
for c in contigs_dict:
	if c in seeds_set:
		contigs_dict[c]['Seed'] = 1
	else:
		contigs_dict[c]['Seed'] = 0

	
	#GC_total += contigs_dict[c]['GC_cont']*contigs_dict[c]['Length']
	#ln_total += contigs_dict[c]['Length']
	
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

contigs_dict = plasmids_preprocessing.get_gene_coverage(mapping_file, contigs_dict)
#contigs_dict['1']['Gene_coverage'] = 0.8333
#contigs_dict['2']['Gene_coverage'] = 0.4163
#contigs_dict['3']['Gene_coverage'] = 0.9163
#contigs_dict['4']['Gene_coverage'] = 0.3333

UBD_rd = 0
UBD_GC = 0
UBD_gd = 0
for c in contigs_dict:
	ln_total += contigs_dict[c]['Length']
	UBD_rd = max(UBD_rd, contigs_dict[c]['Read_depth'])
	UBD_GC = max(UBD_GC, contigs_dict[c]['GC_cont'])
	UBD_gd = max(UBD_gd, contigs_dict[c]['Gene_coverage'])
	print(c, contigs_dict[c]['Gene_coverage'], contigs_dict[c]['Length'])	
print(UBD_rd, UBD_gd, UBD_GC, ln_total)


output_file = open(os.path.join(output_folder, output_filename), "w")
questionable_file = open(os.path.join(output_folder, questionable_filename), "w")
score_file = open(os.path.join(output_folder, score_filename), "w")
contigs_file = open(os.path.join(output_folder, contigs_filename), "w")
logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"w")

n_iter = 0
q_iter = 0

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


	contigs = {}
	contigs_ext = {}
	contigs, contigs_ext = plasmids_preprocessing.contig_vars(m, contigs_dict, contigs, contigs_ext, nplasmids)

	links = {}
	links = plasmids_preprocessing.link_vars(m, links_list, links, nplasmids)

	mean_rd, rd_diff, counted_rd_diff, wtd_rd_diff, counted_wtd_rd_diff = {}, {}, {}, {}, {}	
	mean_rd, rd_diff, counted_rd_diff, wtd_rd_diff, counted_wtd_rd_diff \
		= plasmids_preprocessing.rd_vars(m, contigs_dict, mean_rd, rd_diff, counted_rd_diff, wtd_rd_diff, counted_wtd_rd_diff, nplasmids)

	mean_GC, GC_diff, counted_GC_diff, wtd_GC_diff, counted_wtd_GC_diff = {}, {}, {}, {}, {}	
	mean_GC, GC_diff, counted_GC_diff, wtd_GC_diff, counted_wtd_GC_diff \
		= plasmids_preprocessing.rd_vars(m, contigs_dict, mean_GC, GC_diff, counted_GC_diff, wtd_GC_diff, counted_wtd_GC_diff, nplasmids)

	counted_ln = {}
	counted_ln = plasmids_preprocessing.ln_vars(m, contigs_dict, counted_ln, nplasmids)

	counted_seed = {}
	counted_seed = plasmids_preprocessing.seed_vars(m, contigs_dict, counted_seed, nplasmids)

	wtd_gd, counted_wtd_gd = {}, {}
	wtd_gd, counted_wtd_gd \
		= plasmids_preprocessing.gd_vars(m, contigs_dict, wtd_gd, counted_wtd_gd, nplasmids)

	

	#-----------------------------------------------
	#Setting up the expression for the objective function
	expr = LinExpr()
	alpha1, alpha2, alpha3 = float(alpha1), float(alpha2), float(alpha3) 

	for p in rd_diff:
		#for c in rd_diff[p]:
			#expr.addTerms(alpha1, counted_wtd_rd_diff[p][c])

			#expr.addTerms(alpha2, contigs[p][c])

			#expr.addTerms(-alpha2, counted_wtd_gd[p][c])

			#expr.addTerms(alpha3, counted_wtd_GC_diff[p][c])
		expr.addTerms(alpha1, wtd_rd_diff[p])
		expr.addTerms(-alpha2, wtd_gd[p])
		expr.addTerms(alpha3, wtd_GC_diff[p])	
	#expr.addConstant(2)

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
	#for c in contigs_dict:
	#	expr = LinExpr()
	#	for p in rd:
	#		expr.addTerms(1, rd[p][c])
	#	m.addConstr(expr <= contigs_dict[c]['Read_depth'], "read-depth_cap")
	#	constraint_count += 1		
	
	#Constraint type 6
	#counted-read-depth = read-depth * existence-of-a-contig-in-plasmid
	#for p in counted_rd:
	#	for c in counted_rd[p]:
	#		if c in contigs[p]:
	#			m.addConstr(counted_rd[p][c] <= contigs_dict[c]['Read_depth']*contigs[p][c], "counted_read-depth_1")
	#			m.addConstr(counted_rd[p][c] <= rd[p][c], "counted_read-depth_2")
	#			m.addConstr(counted_rd[p][c] >= rd[p][c] - (1 - contigs[p][c])*contigs_dict[c]['Read_depth'], "counted_read-depth_3")
	#			m.addConstr(counted_rd[p][c] >= 0, "counted_read-depth_4")
	#			constraint_count += 4				
	
	
	#Constraint type 6
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
		m.addConstr(expr <= 1750000, "plasmid_length_upper_bound")
		constraint_count += 1	

	#Constraint type 7
	#counted-seed = seed * existence-of-a-contig-in-plasmid
	for p in counted_seed:
		for c in counted_seed[p]:
			if c in contigs[p]:
				m.addConstr(counted_seed[p][c] == contigs_dict[c]['Seed']*contigs[p][c], "counted_seed")
				constraint_count += 1			
	
	#Constraint type 9
	#Computing mean read depth for each plasmid
	#for p in counted_mean_rd:
	#	expr1 = LinExpr()
	#	expr2 = LinExpr()
	#	for c in counted_mean_rd[p]:
	#		if c in contigs[p]:
	#			m.addConstr(counted_mean_rd[p][c] <= 1000000*contigs[p][c], "counted_mean-rd1")
	#			m.addConstr(counted_mean_rd[p][c] <= mean_rd[p], "counted_mean-rd2")
	#			m.addConstr(counted_mean_rd[p][c] >= mean_rd[p] - (1 - contigs[p][c])*1000000, "counted_mean-rd3")
	#			m.addConstr(counted_mean_rd[p][c] >= 0, "counted_mean-rd4")

	#			expr1.addTerms(contigs_dict[c]['Length'], counted_mean_rd[p][c])
	#			expr2.addTerms(contigs_dict[c]['Length'], counted_rd[p][c])
	#			m.addConstr(expr1 == expr2, "computing-mean-rd")
	#			constraint_count += 5	
	

	#Constraint type 8


	#Constraint type 9
	#Computing weighted read depth diff for each plasmid
	for p in counted_wtd_rd_diff:
		expr1 = LinExpr()
		expr2 = LinExpr()
		for c in counted_wtd_rd_diff[p]:
			expr1.addTerms(contigs_dict[c]['Length'], counted_wtd_rd_diff[p][c])
			expr2.addTerms(contigs_dict[c]['Length'], counted_rd_diff[p][c])
		m.addConstr(expr1 == expr2, "computing-wtd-rd-diff")
		constraint_count += 1

	for p in counted_wtd_GC_diff:
		expr1 = LinExpr()
		expr2 = LinExpr()		
		for c in counted_wtd_GC_diff[p]:
			expr1.addTerms(contigs_dict[c]['Length'], counted_wtd_GC_diff[p][c])
			expr2.addTerms(contigs_dict[c]['Length'], counted_GC_diff[p][c])
		m.addConstr(expr1 == expr2, "computing-wtd-GC-diff")
		constraint_count += 1

	for p in counted_wtd_gd:
		expr1 = LinExpr()
		expr2 = LinExpr()		
		for c in counted_wtd_gd[p]:
			expr1.addTerms(contigs_dict[c]['Length'], counted_wtd_gd[p][c])
			expr2.addTerms(contigs_dict[c]['Gene_coverage'], counted_ln[p][c])
		m.addConstr(expr1 == expr2, "computing-wtd-gd")
		constraint_count += 1			

	#Constraint type 10
	#Absolute value (diff) of mean read depth and plasmid read-depth
	for p in rd_diff:
		for c in rd_diff[p]:
			if c in rd_diff[p]:
				m.addConstr(rd_diff[p][c] >= mean_rd[p] - contigs_dict[c]['Read_depth']*contigs[p][c], "rd-diff_lbd1")
				m.addConstr(rd_diff[p][c] >= contigs_dict[c]['Read_depth']*contigs[p][c] - mean_rd[p], "rd-diff_lbd2")
				constraint_count += 2

	for p in GC_diff:
		for c in GC_diff[p]:
			if c in GC_diff[p]:
				m.addConstr(GC_diff[p][c] >= mean_GC[p] - contigs_dict[c]['GC_cont']*contigs[p][c], "GC-diff_lbd1")
				m.addConstr(GC_diff[p][c] >= contigs_dict[c]['GC_cont']*contigs[p][c] - mean_GC[p], "GC-diff_lbd2")
				constraint_count += 2

	for p in counted_rd_diff:
		for c in counted_rd_diff[p]:
			m.addConstr(counted_rd_diff[p][c] <= UBD_rd*contigs[p][c], "counted-rd-diff1")
			m.addConstr(counted_rd_diff[p][c] <= rd_diff[p][c], "counted-rd-diff2")
			m.addConstr(counted_rd_diff[p][c] >= rd_diff[p][c] - (1 - contigs[p][c])*UBD_rd, "counted-rd-diff3")
			m.addConstr(counted_rd_diff[p][c] >= 0, "counted-rd-diff4")
			constraint_count += 4

	for p in counted_GC_diff:
		for c in counted_GC_diff[p]:
			m.addConstr(counted_GC_diff[p][c] <= UBD_GC*contigs[p][c], "counted-GC-diff1")
			m.addConstr(counted_GC_diff[p][c] <= GC_diff[p][c], "counted-GC-diff2")
			m.addConstr(counted_GC_diff[p][c] >= GC_diff[p][c] - (1 - contigs[p][c])*UBD_GC, "counted-GC-diff3")
			m.addConstr(counted_GC_diff[p][c] >= 0, "counted-GC-diff4")
			constraint_count += 4

	for p in counted_wtd_rd_diff:
		for c in counted_wtd_rd_diff[p]:
			m.addConstr(counted_wtd_rd_diff[p][c] <= UBD_rd*contigs[p][c], "counted-wtd-rd-diff1")
			m.addConstr(counted_wtd_rd_diff[p][c] <= wtd_rd_diff[p], "counted-wtd-rd-diff2")
			m.addConstr(counted_wtd_rd_diff[p][c] >= wtd_rd_diff[p] - (1 - contigs[p][c])*UBD_rd, "counted-wtd-rd-diff3")
			m.addConstr(counted_wtd_rd_diff[p][c] >= 0, "counted-wtd-rd-diff4")
			constraint_count += 4

	for p in counted_wtd_GC_diff:
		for c in counted_wtd_GC_diff[p]:
			m.addConstr(counted_wtd_GC_diff[p][c] <= UBD_GC*contigs[p][c], "counted-wtd-GC-diff1")
			m.addConstr(counted_wtd_GC_diff[p][c] <= wtd_GC_diff[p], "counted-wtd-GC-diff2")
			m.addConstr(counted_wtd_GC_diff[p][c] >= wtd_GC_diff[p] - (1 - contigs[p][c])*UBD_GC, "counted-wtd-GC-diff3")
			m.addConstr(counted_wtd_GC_diff[p][c] >= 0, "counted-wtd-GC-diff4")
			constraint_count += 4

	for p in counted_wtd_gd:
		for c in counted_wtd_gd[p]:
			m.addConstr(counted_wtd_gd[p][c] <= UBD_gd*contigs[p][c], "counted-wtd-gd1")
			m.addConstr(counted_wtd_gd[p][c] <= wtd_gd[p], "counted-wtd-gd2")
			m.addConstr(counted_wtd_gd[p][c] >= wtd_gd[p] - (1 - contigs[p][c])*UBD_gd, "counted-wtd-gd3")
			m.addConstr(counted_wtd_gd[p][c] >= 0, "counted-wtd-gd4")
			constraint_count += 4					
				

	#Constraint type 11
	#Each plasmids should have at least one seed
	for p in counted_seed:
		expr = LinExpr()
		for c in counted_seed[p]:
			expr.addTerms(1, counted_seed[p][c])
		m.addConstr(expr >= 1, "seed_existence")
		constraint_count += 1


	#Constraint type 12
	#Preliminary path constraints
	for p in contigs:
		c_expr = LinExpr()
		for c in contigs[p]:

			c_expr.addTerms(1, contigs[p][c])
		l_expr = LinExpr()
		for e in links[p]:
			l_expr.addTerms(1, links[p][e])
		m.addConstr(c_expr == l_expr + 1, "path_simple")	


	#Constraint type 0
	#Non-negativity constraints
	for p in contigs:
		for c in contigs[p]:
			m.addConstr(contigs[p][c] >= 0, "non-negativity")
		for extr in contigs_ext[p]:
			m.addConstr(contigs_ext[p][extr] >= 0, "non-negativity")
		for e in links[p]:
			m.addConstr(links[p][e] >= 0, "non-negativity")
		
		#for c in rd[p]:
		#	m.addConstr(rd[p][c] >= 0, "non-negativity")
		#	m.addConstr(counted_rd[p][c] >= 0, "non-negativity")
		#	m.addConstr(counted_mean_rd[p][c] >= 0, "non-negativity")
			

		m.addConstr(mean_rd[p] >= 0, "non-negativity")
		m.addConstr(mean_GC[p] >= 0, "non-negativity")
		m.addConstr(wtd_gd[p] >= 0, "non-negativity")
		m.addConstr(wtd_rd_diff[p] >= 0, "non-negativity")
		m.addConstr(wtd_GC_diff[p] >= 0, "non-negativity")

		for c in counted_ln[p]:
			m.addConstr(counted_ln[p][c] >= 0, "non-negativity")
			m.addConstr(counted_seed[p][c] >= 0, "non-negativity")
			m.addConstr(rd_diff[p][c] >= 0, "non-negativity")
			m.addConstr(GC_diff[p][c] >= 0, "non-negativity")
			m.addConstr(counted_wtd_gd[p][c] >= 0, "non-negativity")





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

		#solution_seq[p] = plasmids_postprocessing.get_seq(solution_links[p], soln_ext_dict[p], contigs_dict)			

		reached = {}
		#soln_chr_dict = {}
		circ_seq = {}
		c = 1 
		logfile_3.write("Iteration "+ str(i)+"\n")

		#logfile_4 = open(os.path.join(output_folder,'Final_soln.log'),"w")
		circular = 0	

		circular, circ_plasmid, solution_seq[p] = \
			rmcircular.check_if_circular(reached, solution_seq, circ_seq, c, i, soln_ext_dict, solution_links, links, logfile_3)
		
		if circular == 1:	#If chr is circular, adds corresponding constraint
			x = randint(1, len(circ_seq))
			chosen_seq = set()
			for e in circ_seq[x]['Seq']:
				if e not in chosen_seq:
					chosen_seq.add(e)
			#chosen_seq = circ_seq[x]['Seq']
			print(chosen_seq)
			logfile_3.write("\nChosen seed at iteration "+str(i)+": "+str(x))
			logfile_3.write("\nChosen chr at iteration "+str(i)+"\n")
			expr = LinExpr()
			for e in chosen_seq:
				print("Simple",e)
				if e in links[p]:
					#logfile_3.write(str(species)+": "+ str(adj)+"\n")
					expr.addTerms(1, links[p][e])
				elif e[::-1] in links[p]:
					#logfile_3.write(str(species)+": "+ str(adj[::-1])+"\n")
					expr.addTerms(1, links[p][e[::-1]])
			#logfile_3.write("\n\n")		
			m.addConstr(expr <= len(chosen_seq) - 1, "circular")
						
			print("\n", i, "th circular constraint satisfied\n")
			i += 1

		if i >= 50:
			break	
#_______________

	print('Obj:', m.objVal)

	plasmid_length = {}
	plasmid_gd = {}
	for p in contigs:
		plasmid_length[p] = 0
		plasmid_gd[p] = 0
		for c in contigs[p]:
			if contigs[p][c].x > 0:
				print(contigs[p][c].varName, contigs[p][c].x)
				#print(counted_rd[p][c].varName, counted_rd[p][c].x)
				print(contigs_dict[c]['Read_depth'])
				plasmid_length[p] += contigs_dict[c]['Length']
				plasmid_gd[p] += contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length']
		plasmid_gd[p] = plasmid_gd[p]/plasmid_length[p]
		for e in links[p]:
			if links[p][e].x > 0:
				print(links[p][e].varName, links[p][e].x)		
		print(p, plasmid_length[p])
		print(p, plasmid_gd[p])
		print(mean_rd[p].varName, mean_rd[p].x)
		print("\n")

		if UBD_gd < wtd_gd[p].x:
			print('gd upper bound lower than wtd_gd\n')
			quit()	


	output_file = open(os.path.join(output_folder, output_filename), "a")
	questionable_file = open(os.path.join(output_folder, questionable_filename), "a")
	score_file = open(os.path.join(output_folder, score_filename), "a")
	contigs_file = open(os.path.join(output_folder, contigs_filename), "a")


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
		solution_seq[p], contig_chain = plasmids_postprocessing.get_seq(solution_links[p], soln_ext_dict[p], contigs_dict)
		if plasmid_length[p] >= 1500 and wtd_gd[p].x >= 0.3:
			output_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length[p])+"\t"+"gene_density="+str(plasmid_gd[p])+"\t"+"mean_read_depth="+str(mean_rd[p].x)+"\n")
			output_file.write(solution_seq[p]+"\n")

			contigs_file.write("plasmid_"+str(n_iter)+";")
			for x in contig_chain:
				if x == contig_chain[0]:
					contigs_file.write(x)
				else:
					contigs_file.write(","+x)	
			contigs_file.write("\n")	


			rd_sum, gd_sum, GC_sum = 0, 0, 0
			#for p in rd_diff:
				#for c in rd_diff[p]:
				#	rd_sum += alpha1*counted_wtd_rd_diff[p][c].x
					#expr.addTerms(alpha1*contigs_dict[c]['Length'], counted_diff[p][c])
				#	gd_sum += alpha2*counted_wtd_gd[p][c].x
					#print(c, counted_wtd_gd[p][c].x)
				#	if contigs[p][c].x > 0:
				#		print("Selected: ",c, counted_wtd_gd[p][c].x, counted_wtd_GC_diff[p][c].x, counted_wtd_rd_diff[p][c].x)
					#expr.addTerms(-alpha2*contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length'], counted_rd[p][c])
				#	GC_sum += alpha3*counted_wtd_GC_diff[p][c].x
					#expr.addTerms(-alpha3*(GC_mean-contigs_dict[c]['GC_cont']), counted_ln[p][c])
			#score_file.write("plasmid_"+str(n_iter)+"\t"+str(m.objVal)+"\t"+str(rd_sum)+"\t"+str(gd_sum)+"\t"+str(GC_sum)+"\n")
			score_file.write("plasmid_"+str(n_iter)+"\t"+str(m.objVal)+"\t"+str(wtd_rd_diff[p].x)+"\t"+str(wtd_gd[p].x)+"\t"+str(wtd_GC_diff[p].x)+"\n")
			n_iter += 1
		else:
			questionable_file.write(">plasmid_"+str(q_iter)+"\t"+"length="+str(plasmid_length[p])+"\t"+"gene_density="+str(plasmid_gd[p])+"\t"+"mean_read_depth="+str(mean_rd[p].x)+"\n")	
			questionable_file.write(solution_seq[p]+"\n")
			q_iter += 1	

	#Updating assembly graph and formulation
	for p in contigs:
		for c in contigs[p]:
			if contigs[p][c].x == 1:
				contigs_dict[c]['Read_depth'] = max(0, contigs_dict[c]['Read_depth'] - mean_rd[p].x)

				if c in seeds_set and contigs_dict[c]['Read_depth'] <= 0.5:
					contigs_dict[c]['Seed'] = 0
					seeds_set.remove(c)

				if contigs_dict[c]['Read_depth'] <= 0.05:
					del(contigs_dict[c])				

					for e in links_list:						
						if (c, 'h') in e or (c, 't') in e:
							links_list.remove(e)
			