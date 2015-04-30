#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re

import pandas as pd

import scipy as sp
import scipy.cluster
import numpy as np
from collections import Counter

from file_parse import *
from Rep_sequence_analysis import Clone
from ab_classes import *



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
        Vmut = int(mostCommon(Vreads[cdr3_dict[final_seq_3]].Vmut))
        Jmut = int(mostCommon(Vreads[cdr3_dict[final_seq_3]].Jmut))
        sh = Vmut + Jmut

    return J, final_seq_1, final_seq_2, final_seq_3, num_reads, \
            percent_reads, ABtype, ID_seqs, Vmut, Jmut, sh




def reads_to_clones(args): #reads will need to be self.Reads_slpit_by_V[Vgerm]
                                        #germ_num will be len(self.Reads_split_by_V)
                                        #Vgerm will be the string of the v germ line
    
    all_reads = args[0]
    Vgerm = args[1]
    num_Reads = args[2]
    print args[1]

    reads = all_reads[Vgerm]

    All_clones=[]  
    Clones_split_by_V = {}                    
    Clones_split_by_V[Vgerm] = {}

    #load in cdr3 sequences from Vgerm
    all_cdr3s, cdr3_dict = load_cdrs(reads, 'cdr3')
    all_cdr2s, cdr2_dict = load_cdrs(reads, 'cdr2')
    all_cdr1s, cdr1_dict = load_cdrs(reads, 'cdr1')
    all_IDs, ID_dict = load_cdrs(reads, 'ID')

    #if VDJ didn't find the CDR3, skip to next
    if all_cdr3s == []:
        # Vgerm_complete +=1  
        # print str(Vgerm_complete) +'/' +str(germ_num) + " germlines done!"
        return

    #cluster all cdrs by maximum Hamming length = 1 in each
    T = cluster_into_clones(all_cdr1s, all_cdr2s, all_cdr3s)

    #list of numbers of clones
    clones = [x+1 for x in range(np.amax(T))]

    #find properties of the clone
    for clone_num in clones: #clones is list of numbers from T

        J, final_seq_1, final_seq_2, final_seq_3, \
        num_reads, percent_reads, ABtype, IDs, Vmut, Jmut, sh = find_clone_props(all_cdr1s,\
                                                                all_cdr2s,\
                                                                all_cdr3s, \
                                                                cdr3_dict,\
                                                                T, \
                                                                reads, \
                                                                num_Reads,\
                                                                clone_num,\
                                                                all_IDs)
        #plug it all into a clone object
        All_clones.append(Clone(V=Vgerm, \
                           J=J, \
                           cdr1 = Seq(final_seq_1, generic_protein), \
                           cdr2 = Seq(final_seq_2, generic_protein), \
                           cdr3 = Seq(final_seq_3, generic_protein), \
                           num_reads=num_reads, \
                           percent_reads=percent_reads, \
                           ABtype = ABtype,\
                           IDs = IDs,\
                           Vmut = Vmut,\
                           Jmut = Jmut,\
                           sh = sh))


        Clones_split_by_V[Vgerm][clone_num]=(Clone(V=Vgerm, \
                                            J=J, \
                                            cdr1 = Seq(final_seq_1, generic_protein), \
                                            cdr2 = Seq(final_seq_2, generic_protein), \
                                            cdr3 = Seq(final_seq_3, generic_protein), \
                                            num_reads=num_reads, \
                                            percent_reads=percent_reads,\
                                            ABtype = ABtype,\
                                            IDs = IDs,\
                                            Vmut = Vmut,\
                                            Jmut = Jmut,\
                                            sh = sh))
    
    #print progress 
    # Vgerm_complete +=1  
    # print "%s unique clones found " % np.amax(T)
    # print str(Vgerm_complete) +'/' +str(len(self.Reads_split_by_V)) + " germlines done!\n"

    return All_clones, Clones_split_by_V

def reads_to_clones_star(args):
    # print args[0], args[1], args[2]
    print args[1]
    return reads_to_clones(args[0], args[1], args[2])
