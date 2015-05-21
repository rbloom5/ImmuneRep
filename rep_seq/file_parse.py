#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re
from ab_classes import *
from collections import Counter


def remove_slash(string):
    #removes slashes from names so that they can be used as
    #file names without looking like a directory
    for i in range(len(string)):
        if string[i] == '/':
            string = string[:i] + '%' + string[i+1:]
    return string
            
    

def parse_germ(germ_list):
    germline = []
    mutations = []
    germ_matches = len(germ_list)/3
    
    #pick out the germline name and # of mismatches between sequencer read and germline
    for i in range(germ_matches):
        germline.append(remove_slash(germ_list[i*3]))
        mutations.append(int(germ_list[i*3+2]))
    
    return germline, mutations




def parse_v_j(filepath, num, select_type): #filepath must be file with extension XXX.VDJ.H3.L3.CH1.fa from VDJfasta
	#filename is the .VDJ.H3.L3.CH1.fa to parse
	#number is the max number of reads
	#select_type specifies if you only want one type of Ab (IGHG, IGHM etc.) is optional
	Reads = {}
	num_reads = 0
	if not isinstance(select_type, list):
		select_type = [select_type]
	if not isinstance(filepath, list):
		filepath = [filepath]
	for f in filepath:
		num_reads = 0
		with open(f, "rU") as in_handle:
			v_j=[]

			for line in in_handle:
				if re.search('^>', line):

					#name is the ID in the file
					name = line.split(';')[0].split(' ')[0][1:]
					currentABtype = parse_germ(line.split(';')[6].split(' '))[0]

					#check type
					if select_type and currentABtype and (currentABtype[0] not in select_type):
						continue

					#parse_germ returns the germline and the mismatches between sequencer read and germline
					V, Vmut = parse_germ(line.split(';')[1].split(' '))
					J, Jmut = parse_germ(line.split(';')[3].split(' '))

					#make protein Seq object from cdr3
					cdr3 = Seq(line.split(';')[4], generic_protein)
					cdr2 = Seq(line.split(';')[9], generic_protein)
					cdr1 = Seq(line.split(';')[8], generic_protein)

					#check if all cdrs are present
					if str(cdr1) and str(cdr2) and str(cdr3):
					#fill Reads dict with a Ab_read object with above info
						Reads[name] = Ab_read(name = name, V=V, Vmut=Vmut, J=J, Jmut=Jmut, \
											ABtype = currentABtype, cdr3=cdr3, cdr2=cdr2, cdr1=cdr1)
						num_reads += 1
						if num_reads == num:
							break
	return Reads

def load_cdrs(clones, cdr):
	all_cdrs=[]
	cdr_dict = {}
	for read in clones:
		if cdr == 'cdr3':
			cdr_list = clones[read].cdr3
		elif cdr == 'cdr2':
			cdr_list = clones[read].cdr2
		elif cdr == 'cdr1':
			cdr_list = clones[read].cdr1
		elif cdr == "ID":
			cdr_list = clones[read].name
		elif cdr == "IDs":
			cdr_list = clones[read].IDs
		else:
			return "error in cdr input name"


		if str(cdr_list):
			all_cdrs.append(str(cdr_list)) ###plug this into cluster_into_clones
			cdr_dict[str(cdr_list)] = read #a dict with cdr3 sequence as key
																		#and ID as value

	return all_cdrs, cdr_dict


def mostCommon(items):
    c = Counter(items)
    mode = c.most_common(1)[0][0]
    return mode


def find_best_id(clone, reads):
	cdr3s = []
	for cid in clone.IDs:
		cdr3s.append(str(reads[cid].cdr3))

	final_cdr3 = mostCommon(cdr3s)
	for cid in clone.IDs:
		if str(reads[cid].cdr3)==final_cdr3:
			best_id = cid
			return best_id



def find_and_write(best_id, f, big_file):
	write = 0
	with open(big_file, 'r') as bf:
		for line in bf:

			if write ==1 and re.search('^>', line):
				return

			if line.split(' ')[0][1:] == best_id:
				write=1

			if write==1:
				f.write(line)





