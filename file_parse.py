#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re


class Ab_read:
	def __init__(self, name = '', V = '', Vmut = '', J='', Jmut='', ABtype = '',cdr3=Seq('')):
		self.name = name
		self.V = V
		self.Vmut = Vmut
		self.J = J
		self.Jmut = Jmut
		self.ABtype = ABtype
		self.cdr3 = cdr3



def remove_slash(string):
    #removes slashes from names so that they can be used as
    #file names without looking like a directory
    for i in range(len(string)):
        if string[i] == '/':
            string = string[:i] + '_' + string[i+1:]
    return string
            
    

def parse_germ(germ_list):
    germline = []
    mutations = []
    germ_matches = len(germ_list)/3
    
    #pick out the germline name and # of mismatches between sequencer read and germline
    for i in range(germ_matches):
        germline.append(remove_slash(germ_list[i*3]))
        mutations.append(germ_list[i*3+2])
    
    return germline, mutations




def parse_v_j(filepath, num, select_type): #filepath must be file with extension XXX.VDJ.H3.L3.CH1.fa from VDJfasta
	Reads = {}
	num_reads = 0
	if not isinstance(select_type, list):
		select_type = [select_type]
	with open(filepath, "rU") as in_handle:
		v_j=[]

		for line in in_handle:
			if re.search('^>', line):
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

				
				#fill Reads dict with a Rep_read object with above info
				Reads[name] = Ab_read(name = name, V=V, Vmut=Vmut, J=J, Jmut=Jmut, ABtype = currentABtype, cdr3=cdr3)
				num_reads += 1
				if num_reads == num:
					return Reads
	return Reads

def load_cdr3s(clones):
	all_cdr3s=[]
	cdr3_dict = {}
	for read in clones:
		if str(clones[read].cdr3):
			all_cdr3s.append(str(clones[read].cdr3)) ###plug this into cluster_into_clones
			cdr3_dict[str(clones[read].cdr3)] = read #a dict with cdr3 sequence as key
																		#and ID as value

	return all_cdr3s, cdr3_dict








