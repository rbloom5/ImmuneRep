#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re


class Ab_read:
	def __init__(self, name = '', V = '', Vmut = '', J='', Jmut='', ABtype = '',cdr3=Seq(''),cdr2=Seq(''),cdr1=Seq('')):
		self.name = name
		self.V = V
		self.Vmut = Vmut
		self.J = J
		self.Jmut = Jmut
		self.ABtype = ABtype
		self.cdr3 = cdr3
		self.cdr2 = cdr2
		self.cdr1 = cdr1



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
	if not isinstance(filepath, list):
		filepath = [filepath]
	for f in filepath:
		num_reads = 0
		with open(f, "rU") as in_handle:
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








