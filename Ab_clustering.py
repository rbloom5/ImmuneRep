#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re

import scipy as sp
import scipy.cluster
import numpy as np
from collections import Counter





def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        return len(s1)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def pdist(X,metric):
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            dm[k] = metric(X[i], X[j])
            k += 1
            
            # if k/float(5000000) == k/5000000:
            #     print k/float(len(dm))
    return dm



def cluster_into_clones(seqs,cutoff=1.5,linkage='single'):
    # check trivial cases
    if len(seqs) == 0:
        return np.array([])
        # raise Exception, "chains has nothing it"
    
    #unique_seqs = list(set(seqs))
    seq_idxs = dict( [(j,i) for (i,j) in enumerate(seqs)] )
    
    # check trivial case
    if len(seqs) == 1:
        T = np.array([1]*len(seqs))
        return T
    
    # compute the distance matrix
    Y = pdist(seqs, hamming_distance)
    
    # compute the linkage
    Z = sp.cluster.hierarchy.linkage(Y,method=linkage)
    
    # determine the clusters at level cutoff
    T = sp.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')
    
    return T


def mostCommon(items):
    c = Counter(items)
    mode = c.most_common(1)[0][0]
    return mode


def order_clones(All_clones):
	sorted_clones = {}
	freq_dict={}
	for clone_num in All_clones:
		freq_dict[clone_num] = All_clones[clone_num].num_reads
		sorted_freqs = sorted(freq_dict, key=freq_dict.__getitem__, reverse=True)

		for i in range(len(sorted_freqs)):
			sorted_clones[i+1] = All_clones[sorted_freqs[i]]

	return sorted_clones


def find_clone_props(all_cdr3s, cdr3_dict, T, Vreads, num_Reads, clone_num):
	ABtype = []
	percent_reads = ''
	num_reads = ''
	final_seq = ''
	J=''
	V=''

	read_inds = np.where(T==clone_num)[0] #finds the indices in T which have this clone number
	num_reads = len(read_inds)
	percent_reads = num_reads/float(num_Reads)

	clone_seqs=[]
	for i in read_inds:
		clone_seqs.append(all_cdr3s[i])
	final_seq = mostCommon(clone_seqs) #consensus cdr3 sequence

	#put in J
	if Vreads[cdr3_dict[final_seq]].J:
		J = Vreads[cdr3_dict[final_seq]].J[0]

	#put in type
	if Vreads[cdr3_dict[final_seq]].ABtype:
		ABtype = Vreads[cdr3_dict[final_seq]].ABtype[0]

	return J, final_seq, num_reads, percent_reads, ABtype



