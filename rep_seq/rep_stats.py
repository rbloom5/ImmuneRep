#!/usr/bin/python

def find_distribution(clones,num_chunks, num_clones,num_Reads):
	n = num_chunks 
	clone_fractions = []
	if num_clones > n:
		x=range(num_clones)
		num = float(len(x))/n
		l = [ x [i:i + int(num)] for i in range(0, (n-1)*int(num), int(num))]
		#l is a list of lists - each sub list containing the clone indexes in that .1% chunk
		l.append(x[(n-1)*int(num):])

		for chunk in l:
			sum_reads = 0
			for cl in chunk:
				sum_reads+=clones[cl+1].num_reads

			clone_fractions.append(sum_reads/float(num_Reads))
	else:
		for cl in clones:
			clone_fractions.append(Clones[cl].num_reads/float(num_Reads))

	return clone_fractions



def top_clone_fractions(clones, num_Reads):
	top_10 = 1
	top_100 = 1

	top_clone = clones[1].num_reads/float(num_Reads)

	if clones[10]:
		sum_reads = 0
		for i in range(10):
			sum_reads+=clones[i+1].num_reads

		top_10 = sum_reads/float(num_Reads)

	if clones[100]:
		sum_reads = 0
		for i in range(100):
			sum_reads+=clones[i+1].num_reads

		top_100 = sum_reads/float(num_Reads)

	return top_clone, top_10, top_100


def ABtype_fractions(clones, num_clones, num_Reads):
	AB_counts = {'IGHM': 0, 'IGHG': 0, 'IGHA': 0, 'IGHE': 0, 'IGHD':0}
	AB_unique_counts = {'IGHM': 0, 'IGHG': 0, 'IGHA': 0, 'IGHE': 0, 'IGHD':0}
	for cl in clones:
		if clones[cl].ABtype:
			AB_counts[clones[cl].ABtype]+=clones[cl].num_reads
			AB_unique_counts[clones[cl].ABtype]+=1

	for typ in AB_counts:
		AB_counts[typ] = AB_counts[typ]/float(num_Reads)
		AB_unique_counts[typ] = AB_unique_counts[typ]/float(num_clones)

	return AB_counts, AB_unique_counts




