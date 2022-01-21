from gurobipy import *
import collections
import random
import copy

def other(ext):
	c = 't' if ext == 'h' else 'h'
	return c	

def get_seq(soln_links, neighbour_dict, contigs_dict, contigs):
	plasmid_links = copy.deepcopy(soln_links)
	n_dict = copy.deepcopy(neighbour_dict)

	reached_links = {}
	#reached_contigs = {}
	plasmid = []

	extr_degree = {}
	for p in contigs:
		for c in contigs[p]:
			extr_degree[c] = {}
			if (c,'h') in n_dict:
				extr_degree[c][(c,'h')] = len(n_dict[(c,'h')])
			else:
				extr_degree[c][(c,'h')] = 0
			if (c,'t') in n_dict:		
				extr_degree[c][(c,'t')] = len(n_dict[(c,'t')])	
			else:
				extr_degree[c][(c,'t')] = 0				
	plasmid_ends = []
	for c in extr_degree:
		if extr_degree[c][(c,'h')] > extr_degree[c][(c,'t')]:
			plasmid_ends.append((c,'h'))
		elif extr_degree[c][(c,'h')] < extr_degree[c][(c,'t')]:
			plasmid_ends.append((c,'t'))
	#print(plasmid_ends)
	#print(extr_degree)


	if len(plasmid_ends) == 0:
		linear_found = 1
	else:
		linear_found = 0	

	while len(plasmid_links) != 0:
		plasmid_seg = collections.deque()
		if linear_found == 0:
			#Initial assignment
			curr_vertex = plasmid_ends[0]
			c1 = curr_vertex[0]
			plasmid_seg.append((c1, '+')) if curr_vertex[1] == 'h' else plasmid_seg.append((c1, '-'))
			#plasmid_seg.append((c2, '+')) if neighbour[1] == 't' else plasmid_seg.append((c2, '-'))			
			
			while extr_degree[c1][curr_vertex] != 0: 
				#print(curr_vertex, n_dict[curr_vertex])

				#print(curr_vertex, n_dict)
				neighbour = n_dict[curr_vertex][0]
				n_dict[curr_vertex].remove(neighbour)
				n_dict[neighbour].remove(curr_vertex)

				#print("Neib", neighbour_dict)
				#print("N", n_dict)
				#print(curr_vertex, neighbour, n_dict[curr_vertex])
				if len(n_dict[curr_vertex]) == 0:
					del(n_dict[curr_vertex])
				if len(n_dict[neighbour]) == 0:
					del(n_dict[neighbour])	

				curr_link = (curr_vertex, neighbour)
				if curr_link not in plasmid_links:
					curr_link = curr_link[::-1]
				
				#print(curr_link)
				#print(plasmid_links)
				plasmid_links.remove(curr_link)


				
				ext1, ext2 = curr_link[0], curr_link[1]
				extr_degree[ext1[0]][ext1] -= 1
				extr_degree[ext2[0]][ext2] -= 1

				c1 = neighbour[0]
				plasmid_seg.append((c1, '+')) if neighbour[1] == 't' else plasmid_seg.append((c1, '-'))
				curr_vertex = (c1, other(neighbour[1]))

			#print(plasmid_seg)	
			plasmid.append(plasmid_seg)	
			linear_found = 1	


		else:
			#print(len(n_dict), n_dict)
			curr_vertex = random.choice(list(n_dict))		
			c1 = curr_vertex[0]
			plasmid_seg.append((c1, '+')) if curr_vertex[1] == 'h' else plasmid_seg.append((c1, '-'))
			#plasmid_seg.append((c2, '+')) if neighbour[1] == 't' else plasmid_seg.append((c2, '-'))			
			
			while extr_degree[c1][curr_vertex] != 0: 
				#print(curr_vertex, n_dict[curr_vertex])

				#print(curr_vertex, n_dict)
				neighbour = n_dict[curr_vertex][0]
				n_dict[curr_vertex].remove(neighbour)
				n_dict[neighbour].remove(curr_vertex)
				#print(curr_vertex, neighbour, n_dict[curr_vertex])
				if len(n_dict[curr_vertex]) == 0:
					del(n_dict[curr_vertex])
				if len(n_dict[neighbour]) == 0:
					del(n_dict[neighbour])	

				curr_link = (curr_vertex, neighbour)
				if curr_link not in plasmid_links:
					curr_link = curr_link[::-1]
				
				#print(curr_link)
				#print(plasmid_links)
				plasmid_links.remove(curr_link)


				
				ext1, ext2 = curr_link[0], curr_link[1]
				extr_degree[ext1[0]][ext1] -= 1
				extr_degree[ext2[0]][ext2] -= 1

				c1 = neighbour[0]
				plasmid_seg.append((c1, '+')) if neighbour[1] == 't' else plasmid_seg.append((c1, '-'))
				curr_vertex = (c1, other(neighbour[1]))

			#print(plasmid_seg)	
			plasmid.append(plasmid_seg)
	plasmid_seq = ''
	contig_chain = []						
	for seg in plasmid:
		for contig in seg:
			c, sign = contig[0], contig[1]
			#print("This",c)
			if c in contigs_dict:
				if sign == '+':
					contig_chain.append(str(c)+'+')
					plasmid_seq += contigs_dict[c]['Sequence']
				else:
					plasmid_seq += contigs_dict[c]['Sequence'][::-1]	
					contig_chain.append(str(c)+'-')
	#print("Comparing lists")				
	#print(soln_links)
	#print(plasmid_links)	


	#print("Comparing dicts")				
	#print(neighbour_dict)
	#print(n_dict)					
	return plasmid_seq, contig_chain
		

def trim_zero_end(contig_chain, contigs_dict, end_pt):
	l , r = 0, 0
	if end_pt == 'l':
		c = contig_chain[0][:-1]
		den = contigs_dict[c]['Density']
		if den == 0:
			contig_chain = contig_chain[1:]
			l = 0
		else:
			l = 1
	else:
		c = contig_chain[-1][:-1]
		den = contigs_dict[c]['Density']
		if den == 0:
			contig_chain = contig_chain[:-1]
			r = 0
		else:
			r = 1
	print(contig_chain)						
	return(contig_chain, l, r)			

def refine_chain(contig_chain, contigs_dict):
	l, r = 0, 0
	end_pt = 'l'
	while l == 0:
		contig_chain, l, r = trim_zero_end(contig_chain, contigs_dict, end_pt)
	end_pt = 'r'
	while r == 0:
		contig_chain, l, r = trim_zero_end(contig_chain, contigs_dict, end_pt)	

	length, gd, rd = 0, 0, 0
	seq = ''
	for x in contig_chain:
		c, sign = x[:-1], x[-1]
		length += contigs_dict[c]['Length']
		gd += contigs_dict[c]['Gene_coverage']*contigs_dict[c]['Length']
		if contigs_dict[c]['Read_depth'] < 1:
			rd += contigs_dict[c]['Read_depth']*contigs_dict[c]['Length']
		else:
			rd += 1*contigs_dict[c]['Length']
		if sign == '+':	
			seq += contigs_dict[c]['Sequence']	
		else:
			seq += contigs_dict[c]['Sequence'][::-1]
	gd = gd/length
	rd = rd/length				
	return contig_chain, length, gd, rd, seq