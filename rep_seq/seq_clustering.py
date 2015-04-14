#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re

import pandas as pd

import scipy as sp
import scipy.cluster
import numpy as np
from collections import Counter





def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        return len(s1)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def pdist(X,Y,Z,metric):
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            diff_1 = metric(X[i], X[j])
            diff_2 = metric(Y[i], Y[j])
            diff_3 = metric(Z[i], Z[j])

            dm[k] = max([diff_1, diff_2, diff_3]) #metric(X[i], X[j])
            k += 1
            
            # if k/float(5000000) == k/5000000:
            #     print k/float(len(dm))
    return dm





def cluster_into_clones(cdr_1seqs,cdr_2seqs,cdr_3seqs,cutoff=1.5,linkage='single'):
    # check trivial cases
    if len(cdr_3seqs) == 0:
        return np.array([])
        # raise Exception, "chains has nothing it"
    
    # check trivial case
    if len(cdr_3seqs) == 1:
        T = np.array([1]*len(cdr_3seqs))
        return T
    
    # compute the distance matrix
    #Y = pdist(seqs, hamming_distance)
    Y = pdist(cdr_1seqs, cdr_2seqs, cdr_3seqs, hamming_distance)
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


def find_clone_props(all_cdr1s, all_cdr2s, all_cdr3s, cdr3_dict, T, Vreads, num_Reads, clone_num, all_IDs):
    ABtype = []
    percent_reads = 0
    num_reads = 0
    final_seq_1 = ''
    final_seq_2 = ''
    final_seq_3 = ''
    J=''
    V=''
    Vmut = 0
    Jmut = 0
    sh = 0

    read_inds = np.where(T==clone_num)[0] #finds the indices in T which have this clone number
    num_reads = len(read_inds)
    percent_reads = num_reads/float(num_Reads)

    clone_seqs=[]
    for i in read_inds:
        clone_seqs.append(all_cdr1s[i])
    final_seq_1 = mostCommon(clone_seqs) #consensus cdr1 sequence

    clone_seqs=[]
    for i in read_inds:
        clone_seqs.append(all_cdr2s[i])
    final_seq_2 = mostCommon(clone_seqs) #consensus cdr2 sequence

    clone_seqs=[]
    for i in read_inds:
        clone_seqs.append(all_cdr3s[i])
    final_seq_3 = mostCommon(clone_seqs) #consensus cdr3 sequence

    ID_seqs = []
    for i in read_inds:
        ID_seqs.append(all_IDs[i])

    #put in J
    if Vreads[cdr3_dict[final_seq_3]].J:
        J = Vreads[cdr3_dict[final_seq_3]].J[0]

    #put in type
    if Vreads[cdr3_dict[final_seq_3]].ABtype:
        ABtype = Vreads[cdr3_dict[final_seq_3]].ABtype[0]

    #put in SH
    if Vreads[cdr3_dict[final_seq_3]].Vmut and Vreads[cdr3_dict[final_seq_3]].Jmut:
        Vmut = int(Vreads[cdr3_dict[final_seq_3]].Vmut[0])
        Jmut = int(Vreads[cdr3_dict[final_seq_3]].Jmut[0])
        sh = Vmut + Jmut

    return J, final_seq_1, final_seq_2, final_seq_3, num_reads, \
            percent_reads, ABtype, ID_seqs, Vmut, Jmut, sh







