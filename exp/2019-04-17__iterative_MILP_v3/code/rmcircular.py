from gurobipy import *
import collections

def get_link_soln(soln_ext_dict, solution_links, links):
	for plasmid in links:	#Separating those with pva = 1
		soln_ext_dict[plasmid] = {}
		solution_links[plasmid] = []
		for link in links[plasmid]:
			var = links[plasmid][link]
			if var.x > 0:
				solution_links[plasmid].append(link) 
				soln_ext_dict[plasmid][link[0]] = link[1]
				soln_ext_dict[plasmid][link[1]] = link[0]
	return soln_ext_dict, solution_links

def other(ext):
	c = 't' if ext == 'h' else 'h'
	return c

def get_next_link(contig, ext, ext_dict):
	next_contig, next_ext, next_link = None, None, None
	other_ext = other(ext)[0]
	l = (contig, other_ext)
	if l in ext_dict:
		next_contig = ext_dict[l][0]
		next_ext = ext_dict[l][1]
	if next_contig != None:
		r = (next_contig, next_ext)
		next_link = (l,r)
	return next_link, next_contig, next_ext

def check_if_circular(reached, solution_seq, circ_seq, c, i, soln_ext_dict, solution_links, links, logfile_3):
	circular = 0
	for plasmid in solution_links:
		reached[plasmid] = {}
		solution_seq[plasmid] = {}
		seq_idx = 0

		for link in solution_links[plasmid]:
			#circular = 0
			l, r = link[0], link[1]
			c1, ext1, c2, ext2 = l[0], l[1], r[0], r[1]
			if c1 not in reached[plasmid] and c2 not in reached[plasmid]:		#Every gene would be a part of exactly one chr
				seq_idx += 1
				solution_seq[plasmid][seq_idx] = {}
				d = collections.deque()			#List (queue) of (gene, orientation) tuples
				current_seq = collections.deque()	#List (queue) of adjacecencies

				reached[plasmid][c1] = 1
				reached[plasmid][c2] = 1
		
				d.append((c1, '+')) if ext1 == 'h' else d.append((c1, '-'))
				d.append((c2, '+')) if ext2 == 't' else d.append((c2, '-'))
					
				current_seq.append(link)

				end_reached = 0

				#Follow g2 dir
				contig, ext = c2, ext2
				while end_reached == 0:
					next_link = None
					next_link, next_contig, next_ext = get_next_link(contig, ext, soln_ext_dict[plasmid])

					if next_contig != None and next_contig in reached[plasmid]:
						circular = 1
						seq_type = 'C'
						current_seq.append(next_link)

						circ_seq[c] = {}
						#circ_chr[c]['Species'] = species
						circ_seq[c]['Seq'] = current_seq
						logfile_3.write(str(i)+"\t"+str(current_seq)+"\n")
						c += 1						

						break	#As circular chrom reached. Breaks from while.

					elif next_contig != None:
						if next_ext == 't':
							d.append((next_contig, '+'))
						else:
							d.append((next_contig, '-'))
						current_seq.append(next_link)

						reached[plasmid][next_contig] = 1
						contig = next_contig
						ext = next_ext	

					else:
						end_reached = 1	#As telomere reached

				#Follow g1 dirn if not circular					
				if 	end_reached == 1:
					seq_type = 'L'
					
					contig, ext = c1, ext1
					while end_reached == 1:
						next_link = None
						next_link, next_contig, next_ext = get_next_link(contig, ext, soln_ext_dict[plasmid])

						if next_contig != None:
							if next_ext == 'h':
								d.appendleft((next_contig, '+'))
							else:
								d.appendleft((next_contig, '-'))
							reached[plasmid][next_contig] = 1
							contig = next_contig
							ext = next_ext

						else:
							end_reached = 2

				solution_seq[plasmid][seq_idx]['Seq'] = d
				solution_seq[plasmid][seq_idx]['Type'] = seq_type
				'''
				for adj in d:
					if adj in pva_dict[species]:
						logfile_4.write(str(species)+": "+ str(adj)+"\n")
					elif adj[::-1] in pva_dict[species]:
						logfile_4.write(str(species)+": "+ str(adj[::-1])+"\n")	
				'''		

	#logfile_4.close()
	return circular, circ_seq, solution_seq		