#! /usr/bin/env python

from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import pandas as pd
import numpy as np
import os
import sys
import subprocess
import matplotlib.pyplot as plot

import random
import sys
sys.path.append('/Users/Ryan/SoftwareProjects/MLfunctions')

import cv_averages 
reload(cv_averages)
from cv_averages import cv_metrics

from sklearn.svm import SVC
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier



def getsequence(entry, out_path):
    aaseqs = []
    ntseqs = []
    noHchain = []
    """
    Get the sequence files     
    """
    for e in entry:
        ntSeqObj = []
        RecordExists = False
        translateDict = {"F" : "TTC",
                         "L" : "TTA",
                         "S" : "TCA",
                         "Y" : "TAC",
                         "C" : "TGT",
                         "P" : "CCT",
                         "H" : "CAT",
                         "Q" : "CAG",
                         "R" : "CGC",
                         "I" : "ATC",
                         "M" : "ATG",
                         "T" : "ACA",
                         "N" : "AAC",
                         "K" : "AAG",
                         "V" : "GTT",
                         "A" : "GCA",
                         "D" : "GAC",
                         "E" : "GAA",
                         "G" : "GGC",
                         "W" : "TGG",
                         "X" : "GGT"}
        
        out_file = os.path.join( out_path, "%s_raw.fasta"%e)
        out_file, poo = urllib.urlretrieve("http://opig.stats.ox.ac.uk/webapps/abdb/web_front/data/entries/%s/sequences/%s_raw.fa"%(e,e))
        if os.path.isfile(out_file):       
            with open(out_file, "r") as Retrieved, open("aaABfasta.fasta", 'a+b') as aawritefile, open("ntABfasta.fasta", 'a+b') as ntwritefile:
                for record in SeqIO.parse(Retrieved, "fasta"):
                    if record.description[-1] == 'H':
                        SeqIO.write(record, aawritefile, "fasta")
                        Seqobj = record.seq
                        for i in str(Seqobj):
                            ntSeqObj.append(translateDict[i].lower())
                        ntSeqObj = "".join(ntSeqObj)
                        record.seq = Seq(ntSeqObj) 
                        SeqIO.write(record, ntwritefile, "fasta")                        
                        aaseqs.append(str(Seqobj))
                        ntseqs.append(ntSeqObj)
                        RecordExists = True
                        #print seqs
                        break
            if not RecordExists:
                aaseqs.append(float('NaN'))
                ntseqs.append(float('NaN'))
                noHchain.append(e)
        else:
            aaseqs.append(float('NaN'))
            ntseqs.append(float('NaN'))
            noHchain.append(e)
        os.remove(out_file)
    return {'aa': aaseqs, 'nt': ntseqs, 'noHchain': noHchain}
  


def BlastCompare(seq1, seq2):
	out = []
	#optimized for short sequences, outfmt = XML
	blastOptions = "-evalue=200000 -word_size=2 -matrix='PAM30' -comp_based_stats='0' -outfmt=5"
	query = "-query <(echo -e '>Name\n" + seq1 +"') "
	subject = "-subject <(echo -e '>Name\n" + seq2 +"') "

	blastString = "blastp " + query + subject + blastOptions


	# Run BLAST and parse the output as XML
	process = subprocess.Popen(
		args=blastString,
		stdout=PIPE,
		stderr = subprocess.STDOUT,
		shell=True,
		executable='/bin/bash',
		close_fds=True)
	output=process.communicate()[0]
	blast_result_record = NCBIXML.read(StringIO(output))
    

	#If blast finds a 'hit' record scores, if not, record zeros
	if len(blast_result_record.alignments)>0 :
		for alignment in blast_result_record.alignments:
			for hsp in alignment.hsps[0:1]:
				#save data 
				out.append(alignment.length)
				out.append(hsp.score)
				out.append(hsp.expect)
				out.append(hsp.bits)
				out.append(hsp.align_length)
	else:
		out.append(0)
		out.append(0)
		out.append(0.5) # for evalue - essentially saying the probability the sequences are different is 50-50
		out.append(0)
		out.append(0)





with open('/Users/Ryan/Downloads/summary.tsv', 'rb') as summaryfile:
    pdb = pd.read_csv(summaryfile, sep='\t', header = 0, index_col = 0)

pdb.dropna(subset = ['antigen_name'], inplace = True)
pdb["index"] = pdb.index
pdb.drop_duplicates(cols='index', inplace=True)
del pdb["index"]
indices = pdb.index

#put all the antibodies that bind the same antigen in the same group
grouped = pdb.groupby(by=['antigen_name'])

groupDict = {}
lengths = []
for name, group in grouped:
    groupDict[name] = len(group)
    lengths.append(len(group))

# plt.hist(lengths, bins = np.amax(lengths))
# plt.xlabel('number of antibodies in group')
# plt.ylabel('number of groups')
# plt.show()

s = getsequence(list(indices),'/Users/Ryan/ABsequences/' )






with open('ABsequences/ntABfasta.H3.acc.txt') as f:
    pdb['CDR-H3'] = float('NaN')
    for line in f:
        tempLine = line.split("\t")
        index = tempLine[0][:4]
        pdb['CDR-H3'][index] = tempLine[1][:-1]

pdb.dropna(subset = ['CDR-H3'], inplace = True)





#pairwiseAlign(df, feature)
feature = "CDR-H3"
df = pdb
subsetDiff = 100


#Assign antigen name as index to each antibody to easily check to see if same or different
grouped = df.groupby(by=['antigen_name'])
groupDict = grouped.groups
numGroups = len(groupDict)
#construct new inverted dict
iGD = {}
for index in groupDict:
    for item in groupDict[index]:
        iGD[item] = index

#The features from blast which can be used to compute similarity
outData = {feature+" label" : [],
           "class" : [],
           feature+" score" : [],
           feature+" evalue" : [],
           feature+" bit_score" : [],
           feature+" length" : [], 
           feature+" hsp-length" : []}



count = 0
subsetCount = subsetDiff

#Loop through indices, make pairwise comparisons between sequences
for index in df.index:
    count+=1
    print "progress: " + str( (((len(df.index))**2/2) - ((len(df.index) - count)**2/2)) /float((len(df.index))**2/2)*100-1) + "%"
    for index2 in df.index[count:]:
        
        #if antigen is the same, compare them, if different, compare them if you are at subsetDiff
        if iGD[index] == iGD[index2] or subsetCount == subsetDiff:
            if subsetCount == subsetDiff:
                subsetCount = 1

            #make a label from antibody names    
            outData[feature+" label"].append(index+'-'+index2)

            #assign class to comparison (1 if same, 0 if different)
            if iGD[index] == iGD[index2]:
                outData["class"].append(1)
            else:
                outData["class"].append(0)

            #Compare the two sequences using blastp
            blastScores = BlastCompare(df[feature][index],df[feature][index2])
            
            outData[feature+" length"].append(blastScores[0])
            outData[feature+" score"].append(blastScores[1])
            outData[feature+" evalue"].append(blastScores[2])
            outData[feature+" bit_score"].append(blastScores[3])
            outData[feature+" hsp-length"].append(blastScores[4])

        else:
           subsetCount+=1 

        
print pd.DataFrame(outData)







