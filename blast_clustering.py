#!/usr/bin/python
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
import re

import pandas as pd

import subprocess
from subprocess import PIPE
import scipy as sp
import scipy.cluster
import numpy as np
from collections import Counter
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier

from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment



def blast_score(query_cdrs, subject_cdrs):
    blastOptions = "-evalue=200000 -word_size=2 -matrix='PAM30' -comp_based_stats='0' -outfmt=5"
    outData ={}
    for i in range(3):

        query = "-query <(echo -e '>Name\n" + query_cdrs[i] +"') "
        subject = "-subject <(echo -e '>Name\n" + subject_cdrs[i] +"') "       
        blastString = "blastp " + query + subject + blastOptions


#             # Run BLAST and parse the output as XML
        process = subprocess.Popen(
            args=blastString,
            stdout=PIPE,
            stderr = subprocess.STDOUT,
            shell=True,
            executable='/bin/bash',
            close_fds=True)

        output=process.communicate()[0]
        blast_result_record = NCBIXML.read(StringIO(output))
            


        if len(blast_result_record.alignments)>0 :
            for alignment in blast_result_record.alignments:
                for hsp in alignment.hsps[0:1]:
                    #save data 
                    outData[i] = np.array([hsp.score, hsp.expect, hsp.align_length, alignment.length, hsp.bits])

        else:
            outData[i] = np.array([0, .5, 0, 0, 0])

    return np.concatenate((outData[0], outData[1], outData[2]), axis=1)



def Bio_align(query_cdrs, subject_cdrs):
	outData ={}
	matrix = matlist.pam30
	for i in range(3):
		a = pairwise2.align.globalds(query_cdrs[i], subject_cdrs[i], matrix, -23, -5)[0]
		outData[i] = np.array([a[2], a[4]-a[3]])

	return np.concatenate((outData[0], outData[1], outData[2]), axis=1)






def blast_dist(X):
    rf = joblib.load('/Users/Ryan/SoftwareProjects/ImmuneRep/clustering_model/ABcluster_rf_bio.pkl')
    scaler = joblib.load('/Users/Ryan/SoftwareProjects/ImmuneRep/clustering_model/scaler.pkl')
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            blast_scores = Bio_align(X[i], X[j])
            blast_scores = scaler.transform(blast_scores)	
            dm[k] = rf.predict_proba(blast_scores)[0][0] #metric(X[i], X[j])
            k += 1
            
            # if k/float(5000000) == k/5000000:
            #     print k/float(len(dm))
    return dm






def similarity_cluster(cdrs,proba_cutoff=.5,linkage='average'):
    # check trivial cases
    if len(cdrs) == 0:
        return np.array([])
        # raise Exception, "chains has nothing it"
    
    # check trivial case
    if len(cdrs) == 1:
        T = np.array([1]*len(cdrs))
        return T
    
    # compute the distance matrix
    #Y = pdist(seqs, hamming_distance)
    print "calculating blast scores..."
    Y = blast_dist(cdrs)
    # compute the linkage
    print "computing linkages..."
    Z = sp.cluster.hierarchy.linkage(Y,method=linkage)
    
    # determine the clusters at level cutoff
    print "clustering..."
    T = sp.cluster.hierarchy.fcluster(Z,proba_cutoff,criterion='distance')
    
    return T




def find_cluster_props(all_cdr3s, cdr3_dict, Vreads, read_inds):
    Js=[]
    ABtypes = []
    percent_reads = 0
    num_reads = 0


    #num_reads = len(read_inds)
    #percent_reads = num_reads/float(num_Reads)

    for read_ind in read_inds:

        #put in J
        if Vreads[read_ind+1].J:
            Js.append(Vreads[read_ind+1].J)

        #put in type
        if Vreads[read_ind+1].ABtype:
            ABtypes.append(Vreads[read_ind+1].ABtype)


        percent_reads+=Vreads[read_ind+1].percent_reads
        num_reads+=Vreads[read_ind+1].num_reads



    return Js, ABtypes, num_reads, percent_reads




