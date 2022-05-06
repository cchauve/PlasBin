__author__ = 'amane'

#-------------------

#USAGE: 
#time python plasmids_iterative.py --ag assembly.gfa --map mapping.csv --seeds seed_contigs.csv \
#				 --out output_dir --alpha1 alpha_1 --alpha2 alpha_2 

from gurobipy import *
from sys import argv
import os
import time
from random import randint
import math
import argparse
import plasmids_preprocessing
import plasmids_postprocessing
import rmcircular

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--ag", help="Path to assembly graph file")
	parser.add_argument("--map", help="Path to gene to contig mapping file")
	parser.add_argument("--seeds", help="Path to seed contigs file")
	parser.add_argument("--out", help="Path to output dir")
	parser.add_argument("--alpha1", nargs='?', const = 1, type=int, default = 1, help="Weight of gene density term")
	parser.add_argument("--alpha2", nargs='?', const = 1, type=int, default = 1, help="Weight of GC content term")
	args = parser.parse_args()

	output_dir = args.out
	assembly_file = args.ag
	mapping_file = args.map
	seeds_file = args.seeds
	alpha1 = args.alpha1
	alpha2 = args.alpha2

	#Naming and creating output files
	ratios = str(alpha1) + '.' + str(alpha2)
	output_folder = output_dir + '/' + ratios + '/MILP'
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)
	
	output_fasta = 'putative_plasmids.fasta'
	ques_fasta = 'questionable_plasmids.fasta'
	score_filename = 'MILP_objective.csv'
	output_contigs = 'contig_chains.csv'
	ques_contigs = 'questionable_contig_chains.csv'
	components = 'components.csv'

	#-----------------------------------------------
	#Main program
	contigs_dict = {}
	links_list = []
	seeds_set = set()
	rd_dict = {}
	edge_fam_dict = {}
	contigs_dict, links_list, rd_dict, edge_fam_dict = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list, rd_dict, edge_fam_dict)
	seeds_set = plasmids_preprocessing.get_seeds(seeds_file, seeds_set, rd_dict)
	components_dict = {}
	used_contigs = {}
	#ques_component_dict = {}
	
	ln_total = 0
	n_contigs = 0
	for c in contigs_dict:
		if c in seeds_set:
			contigs_dict[c]['Seed'] = 1
		else:
			contigs_dict[c]['Seed'] = 0
		n_contigs += 1
		ln_total += contigs_dict[c]['Length']		

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

	UBD_rd = 0
	UBD_GC = 0
	UBD_gd = 0
	for c in contigs_dict:
		UBD_rd = max(UBD_rd, contigs_dict[c]['Read_depth'])
		UBD_GC = max(UBD_GC, contigs_dict[c]['GC_cont'])
		UBD_gd = max(UBD_gd, contigs_dict[c]['Gene_coverage'])
	#print(UBD_rd, UBD_gd, UBD_GC, ln_total)

	input_details = 'details.csv'
	details_file = open(os.path.join(output_folder, input_details), "w")
	details_file.write("Contig"+"\t"+'Gene_coverage'+ "\t"+'Read_depth'+ "\t"+'GC_cont'+"\t"+ "Length"+"\n")
	for c in contigs_dict:
		details_file.write(c+"\t"+ str(contigs_dict[c]['Gene_coverage'])+ "\t"+str(contigs_dict[c]['Read_depth'])+ "\t"+str(contigs_dict[c]['GC_cont'])+"\t"+ str(contigs_dict[c]['Length'])+"\n")
		details_file.write("Link list:"+"\t")
		for link in links_list:
			if c == link[0][0] or c == link[1][0]:
				details_file.write(str(link)+"\t")
		details_file.write("\n")
				 
	output_fasta_file = open(os.path.join(output_folder, output_fasta), "w")
	ques_fasta_file = open(os.path.join(output_folder, ques_fasta), "w")
	score_file = open(os.path.join(output_folder, score_filename), "w")
	contigs_file = open(os.path.join(output_folder, output_contigs), "w")
	ques_contigs_file = open(os.path.join(output_folder, ques_contigs), "w")
	components_file = open(os.path.join(output_folder, components), "w")

	logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"w")

	n_iter = 0
	q_iter = 0
	n_comp = 0
	q_comp = 0

	
	while len(seeds_set) > 0:
		
		print("\n\n\n\n\n")
		#-----------------------------------------------
		#Initializing the ILP
		m = Model("Plasmids")
		m.params.LogFile= os.path.join(output_folder,'m.log')
		m.setParam(GRB.Param.TimeLimit, 2400.0)
		m.setParam(GRB.Param.MIPGap, 0.05)

		contigs = {}
		degree = {}
		contigs_ext = {}
		ceil = {}
		contigs, degree, contigs_ext, ceil = plasmids_preprocessing.contig_vars(m, contigs_dict, contigs, degree, contigs_ext, ceil)

		links = {}
		links = plasmids_preprocessing.link_vars(m, links_list, links)

		#rd_diff, counted_rd_diff = {}, {}	
		#rd_diff, counted_rd_diff \
		#	= plasmids_preprocessing.rd_vars(m, contigs_dict, rd_diff, counted_rd_diff)

		mean_GC, counted_GC_mean, GC_diff, counted_GC_diff = {}, {}, {}, {}
		mean_GC, counted_GC_mean, GC_diff, counted_GC_diff \
			= plasmids_preprocessing.GC_vars(m, contigs_dict, mean_GC, counted_GC_mean, GC_diff, counted_GC_diff)

		counted_ln = {}
		counted_ln = plasmids_preprocessing.ln_vars(m, contigs_dict, counted_ln)

		counted_seed = {}
		counted_seed = plasmids_preprocessing.seed_vars(m, contigs_dict, counted_seed)

		#-----------------------------------------------
		#Setting up the expression for the objective function
		expr = LinExpr()
		alpha1, alpha2 = float(alpha1), float(alpha2)

		for p in GC_diff:
			for c in GC_diff[p]:
				expr.addTerms(-alpha1*contigs_dict[c]['Gene_coverage'], contigs[p][c])
				expr.addTerms(alpha2, counted_GC_diff[p][c])
		m.setObjective(expr, GRB.MINIMIZE)

		#-----------------------------------------------
		#Setting up constraints

		constraint_count = 0
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
		#for p in rd_diff:
		#	for c in rd_diff[p]:
		#		if c in rd_diff[p]:
		#			m.addConstr(rd_diff[p][c] >= contigs[p][c] - contigs_dict[c]['Read_depth']/rd_graph, "rd-diff_lbd1")
		#			m.addConstr(rd_diff[p][c] >= 0, "rd-diff_lbd2")
		#			constraint_count += 2
		#for p in counted_rd_diff:
		#	for c in counted_rd_diff[p]:
		#		expr = LinExpr()
		#		expr.addTerms(0.01, contigs[p][c])
		#		m.addConstr(ceil[p][c] <= expr + 0.99, "get-ceil-ubd")
		#		m.addConstr(ceil[p][c] >= expr, "get-ceil-lbd")

		#		m.addConstr(counted_rd_diff[p][c] <= UBD_rd*ceil[p][c], "counted-rd-diff1")
		#		m.addConstr(counted_rd_diff[p][c] <= rd_diff[p][c], "counted-rd-diff2")
		#		m.addConstr(counted_rd_diff[p][c] >= rd_diff[p][c] - (1 - ceil[p][c])*UBD_rd, "counted-rd-diff3")
		#		m.addConstr(counted_rd_diff[p][c] >= 0, "counted-rd-diff4")
		#		constraint_count += 4

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
				#m.addConstr(rd_diff[p][c] >= 0, "non-negativity")
				m.addConstr(GC_diff[p][c] >= 0, "non-negativity")

		#Checking for existence of disconnected cyclic components	
		#Adding cycle elimination constraints if needed	
		logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"a")
		circular = 1
		i = 1
		while circular == 1:
			start = time.time()
			m.optimize()
			stop = time.time()
			duration = stop - start

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

			circular = 0	

			circular, circ_seq, solution_seq[p] = \
				rmcircular.check_if_circular(solution_seq, circ_seq, k, i, soln_ext_dict, solution_links, contigs, logfile_3)
			
			if circular == 1:
				for x in circ_seq:
					chosen_seq = set()	
					chosen_seq = circ_seq[x]['Seq']
				
					expr = LinExpr()
					for e in chosen_seq:
						if e in links[p]:
							expr.addTerms(1, links[p][e])
						elif e[::-1] in links[p]:
							expr.addTerms(1, links[p][e[::-1]])
					logfile_3.write("\n\n")		
					m.addConstr(expr <= len(chosen_seq) - 1, "circular")
							
					i += 1
			if i >= 50:
				break	
	#_______________
	#
		plasmid_length = {}
		plasmid_gd = {}
		plasmid_rd = {}
		for p in contigs:
			plasmid_length[p] = 0
			plasmid_gd[p] = 0
			plasmid_rd[p] = 0
			for c in contigs[p]:
				if contigs[p][c].x > 0:
					plasmid_length[p] += contigs_dict[c]['Length']
					plasmid_gd[p] += contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length']
					#if contigs_dict[c]['Read_depth'] <= rd_graph:
					plasmid_rd[p] += contigs_dict[c]['Read_depth']*contigs_dict[c]['Length']
					#else:
					#	plasmid_rd[p] += rd_graph*contigs_dict[c]['Length']
			plasmid_gd[p] = plasmid_gd[p]/plasmid_length[p]
			plasmid_rd[p] = plasmid_rd[p]/plasmid_length[p]
			#for e in links[p]:
			#	if links[p][e].x > 0:
			#		print(links[p][e].varName, links[p][e].x)
	
		output_fasta_file = open(os.path.join(output_folder, output_fasta), "a")
		ques_fasta_file = open(os.path.join(output_folder, ques_fasta), "a")
		score_file = open(os.path.join(output_folder, score_filename), "a")
		contigs_file = open(os.path.join(output_folder, output_contigs), "a")
		ques_contigs_file = open(os.path.join(output_folder, ques_contigs), "a")
		components_file = open(os.path.join(output_folder, components), "a")

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
						
			contig_chain, plasmid_length[p], plasmid_gd[p], plasmid_rd[p], solution_seq[p] = plasmids_postprocessing.refine_chain(contig_chain, contigs_dict)	

			#Sorting into putative and questionable
			if plasmid_length[p] >= 1500 and plasmid_gd[p] >= 0.3:
				output_fasta_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length[p])+"\t"+"gene_density="+str(plasmid_gd[p])+"\t"+"mean_read_depth="+str(plasmid_rd[p])+"\n")
				if len(solution_seq[p]) != 0:							#if length of contig chain is more than zero
					output_fasta_file.write(solution_seq[p]+"\n")
				else:													#if plasmid is a single contig
					c = contig_chain[0][:-1]								
					output_fasta_file.write(contigs_dict[c]['Sequence']+"\n")
				
				#Printing putative contig chains to file
				contigs_file.write("plasmid_"+str(n_iter)+";")
				components_file.write("plasmid_"+str(n_iter)+";\n")
				for x in contig_chain:
					if x == contig_chain[0]:
						contigs_file.write(x)
					else:
						contigs_file.write(","+x)					
				contigs_file.write("\n")
				if len(contig_chain) == 1:
					components_file.write(contig_chain[0])
					components_file.write("\n")
				#Printing putative contig chain details to file	
				contigs_file.write("#Contig\tGene_coverage\tGC_cont\tRead_depth\tLength\n")		
				for x in contig_chain:
					c = x[:-1]
					contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['Gene_coverage'])+"\t"+str(contigs_dict[c]['GC_cont'])+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
				contigs_file.write("\n")

				for x in plasmid_parts:
					components_dict[n_comp] = x
					n_comp += 1
					for i in range(len(x)):
						y = x[i]
						if i == 0:
							components_file.write(str(y[0]+y[1]))
						else:
							components_file.write(","+str(y[0]+y[1]))	
					components_file.write("\n")

				#Recording scores
				gd_sum, GC_sum = 0, 0
				for p in contigs:
					for c in contigs[p]:
						gd_sum += alpha1*(-contigs_dict[c]['Gene_coverage'])*contigs[p][c].x
						GC_sum += alpha2*counted_GC_diff[p][c].x
				score_file.write("putative_plasmid_"+str(n_iter)+"\t\t"+str(m.objVal)+"\t"+str(gd_sum)+"\t"+str(GC_sum)+"\n")
				n_iter += 1

			else:
				#Same for questionable plasmids		
				ques_fasta_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length[p])+"\t"+"gene_density="+str(plasmid_gd[p])+"\t"+"mean_read_depth="+str(plasmid_rd[p])+"\n")
				if len(solution_seq[p]) != 0:
					ques_fasta_file.write(solution_seq[p]+"\n")
				else:
					c = contig_chain[0][:-1]
					ques_fasta_file.write(contigs_dict[c]['Sequence']+"\n")
		
				ques_contigs_file.write("plasmid_"+str(q_iter)+";")
				for x in contig_chain:
					if x == contig_chain[0]:
						ques_contigs_file.write(x)

					else:
						ques_contigs_file.write(","+x)
				ques_contigs_file.write("\n")		
				for x in contig_chain:
					c = x[:-1]
					ques_contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['Gene_coverage'])+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
				ques_contigs_file.write("\n")

				gd_sum, GC_sum = 0, 0
				for p in contigs:
					for c in contigs[p]:
						gd_sum += alpha1*(-contigs_dict[c]['Gene_coverage'])*contigs[p][c].x
						GC_sum += alpha2*counted_GC_diff[p][c].x
				score_file.write("questionable_plasmid_"+str(n_iter)+"\t"+str(m.objVal)+"\t"+str(gd_sum)+"\t"+str(GC_sum)+"\n")
				q_iter += 1

		#Updating assembly graph and formulation
		for p in contigs:
			for c in contigs[p]:
				if contigs[p][c].x > 0:
					contigs_dict[c]['Read_depth'] = max(0, contigs_dict[c]['Read_depth'] - contigs[p][c].x)
					if c in seeds_set and contigs_dict[c]['Read_depth'] <= 0.5:
						contigs_dict[c]['Seed'] = 0
						seeds_set.remove(c)

					if contigs_dict[c]['Read_depth'] <= 0.05:
						orig_c = c.split('_')[0]
						used_contigs[orig_c] = {}
						used_contigs[orig_c]['Gene_coverage'] = contigs_dict[c]['Gene_coverage']
						used_contigs[orig_c]['Length'] = contigs_dict[c]['Length']
						used_contigs[orig_c]['GC_cont'] = contigs_dict[c]['GC_cont']
						used_contigs[orig_c]['Sequence'] = contigs_dict[c]['Sequence']
						del(contigs_dict[c])				

						temp_list = []
						for e in links_list:						
							if (c, 'h') not in e and (c, 't') not in e:
								temp_list.append(e)
						links_list = temp_list
	#_______________
	#
	#Retaining unique plasmids from output		
	
	component_contigs_list = []
	count = 0
	for k in components_dict: #For every individual plasmid component
		component = components_dict[k] 
		component_contigs = set([x[0].split('_')[0] for x in component])
		#print(component_contigs)
		temp_list = []
		flag = 0 #flag = 1 if set of contigs of current plasmid is subset of set of contigs of any previous plasmids
		for contig_set in component_contigs_list:
			if contig_set[0].issubset(component_contigs) == False: 
				temp_list.append(contig_set)
			elif contig_set[0] == component_contigs:
				temp_list.append(contig_set)
			if flag == 0 and component_contigs.issubset(contig_set[0]):
				flag = 1
		if flag == 0:
			temp_list.append((component_contigs,count))
		component_contigs_list = temp_list	
		count += 1	
		
	for pair in component_contigs_list:
		idx = pair[1]
		temp_len = 0
		temp_seq = []
		gene_free_len = []
		gene_free_seq = []
		chain = component_contigs_list[idx][0]
			
		for contig in chain:

			gd = used_contigs[contig]['Gene_coverage']
			ln = used_contigs[contig]['Length']

			if gd >= 0.3:
				if temp_len > 10000:
					gene_free_len.append(temp_len)
					gene_free_seq.append(temp_seq)
				temp_len = 0
				temp_seq = []
			else:
				temp_len += ln
				temp_seq.append(contig)

		for x in range(len(gene_free_seq)):
			components_list[idx] = plasmids_postprocessing.remove_sub_list(gene_free_seq,chain)			
			
	#Printing unique plasmids to output file
	output_fasta_file = open(os.path.join(output_folder, output_fasta), "w")
	ques_fasta_file = open(os.path.join(output_folder, ques_fasta), "w")
	components_file = open(os.path.join(output_folder, components), "w")
	count = 1
	for pair in component_contigs_list:
		idx = pair[1]
		chain = list(component_contigs_list[idx][0])
		if len(chain) > 1 and chain[0] == chain[-1]:
			chain = chain[:-1]
		length = 0
		gd = 0
		rd = 1
		seq = ''
		components_file.write("plasmid_"+str(count)+";")
		for contig in chain:
			length += used_contigs[contig]['Length']
			gd += used_contigs[contig]['Gene_coverage']*used_contigs[contig]['Length']
			seq += used_contigs[contig]['Sequence']
			components_file.write(str(contig)+",")
		gd = gd/length
		if gd >= 0.3 and length >= 1500: 	
			components_file.write("\n")

			output_fasta_file.write(">plasmid_"+str(count)+"\t"+"length="+str(length)+"\t"+"gene_density="+str(gd)+"\t"+"mean_read_depth="+str(1)+"\n")
			output_fasta_file.write(seq+"\n")
		else:
			ques_fasta_file.write(">plasmid_"+str(count)+"\t"+"length="+str(length)+"\t"+"gene_density="+str(gd)+"\t"+"mean_read_depth="+str(1)+"\n")
			ques_fasta_file.write(seq+"\n")			
		count+=1	
		
		
							


