from gurobipy import *
import collections

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

def get_seq(plasmid_links, ext_dict, contigs_dict):
	reached_links = {}
	reached_contigs = {}
	plasmid = []
	for link in plasmid_links:
		#print(link)
		if link not in reached_links:

			reached_links[link] = 1		
		
			#circular = 0
			l, r = link[0], link[1]
			c1, c2 , ext1, ext2 = l[0], r[0], l[1], r[1]

			if c1 not in reached_contigs and c2 not in reached_contigs:
				plasmid_seg = collections.deque()
				curr_seg = collections.deque()

				reached_contigs[c1] = 1
				reached_contigs[c2] = 1

				plasmid_seg.append((c1, '+')) if ext1 == 'h' else plasmid_seg.append((c1, '-'))
				plasmid_seg.append((c2, '+')) if ext2 == 't' else plasmid_seg.append((c2, '-'))
					
				curr_seg.append(link)

				end_reached = 0

				#Follow c2 direction first
				contig, ext = c2, ext2
				while end_reached == 0:
					next_link = None
					next_link, next_contig, next_ext = get_next_link(contig, ext, ext_dict)

					if next_contig != None and next_contig in reached_contigs:
						#circular = 1
						#seg_type = 'C'
						curr_seg.append(next_link)
						plasmid.append(plasmid_seg)					
						break	#As circular chrom reached. Breaks from while.

					elif next_contig != None:
						if next_ext == 't':
							plasmid_seg.append((next_contig, '+'))
						else:
							plasmid_seg.append((next_contig, '-'))
						curr_seg.append(next_link)
						reached_links[next_link] = 1
						reached_contigs[next_contig] = 1
						contig = next_contig
						ext = next_ext	

					else:
						end_reached = 1	#As telomere reached	

				#Follow c1 dirn if not circular					
				if end_reached == 1:
					#seg_type = 'L'	
					contig, ext = c1, ext1
					while end_reached == 1:
						next_link = None
						next_link, next_contig, next_ext = get_next_link(contig, ext, ext_dict)

						if next_contig != None:
							if next_ext == 'h':
								plasmid_seg.appendleft((next_contig, '+'))
							else:
								plasmid_seg.appendleft((next_contig, '-'))
							reached_links[next_link] = 1	
							reached_contigs[next_contig] = 1
							contig = next_contig
							ext = next_ext

						else:
							plasmid.append(plasmid_seg)
							end_reached = 2

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
	print(contig_chain)				
	return plasmid_seq, contig_chain	 

			
