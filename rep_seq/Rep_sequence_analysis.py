#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re

import numpy as np
from collections import Counter
import itertools
from joblib import Parallel, delayed  
import multiprocessing

import vj_split
reload(vj_split)
from vj_split import *

import file_parse
reload(file_parse)
from file_parse import *

import seq_clustering
reload(seq_clustering)
from seq_clustering import *

import blast_clustering
reload(blast_clustering)
from blast_clustering import *

import rep_stats
reload(rep_stats)
from rep_stats import *


import matplotlib.pyplot as plt
try:
	import seaborn as sns
except:
	pass






class Clone:
    def __init__(self, V = '', J='', cdr3=Seq(''), cdr2=Seq(''), cdr1=Seq(''), \
    			ABtype = '', num_reads = None, percent_reads = None, IDs = [],\
    			Vmut = None, Jmut = None, sh = None):
        self.V = V
        self.J = J
        self.cdr3 = cdr3
        self.cdr2 = cdr2
        self.cdr1 = cdr1
        self.num_reads = num_reads
        self.percent_reads = percent_reads
        self.ABtype = ABtype
        self.IDs = IDs
        self.Vmut = Vmut
        self.Jmut = Jmut
        self.sh = sh

	def __str__(self):
		return u'Clone V={V} J={J}, {num_reads} reads'.format(
			V=self.V,
			J=self.J,
			num_reads=self.num_reads,
		)

	def __repr__(self):
		return self.__str__()




class Cluster:
    def __init__(self, V = '', Js=[], cdr3s=[], cdr2s=[], cdr1s=[], \
    			ABtypes = [], num_reads = None, percent_reads = None, \
    			IDs = []):
        self.V = V
        self.Js = Js
        self.ABtypes = ABtypes
        self.cdr3s = cdr3s
        self.cdr2s = cdr2s
        self.cdr1s = cdr1s
        self.num_reads = num_reads
        self.percent_reads = percent_reads
        self.IDs = IDs

	def __str__(self):
		return u'Cluster V={V}, {num_reads} total reads'.format(
		V=self.V,
		num_reads=self.num_reads,
	)

	def __repr__(self):
		return self.__str__()
 


        

class Rep_seq:

	def __init__(self, filepath, num = False, ABtype = []): #filepath must be file with extension XXX.VDJ.H3.L3.CH1.fa from VDJfasta
		#parses the data from XXX.VDJ.H3.L3.CH1.fa into usable format
		#num is an optional input where you can set the max number of reads to take
		#from the file if you don't want them all

		#self.Reads is a dict of Ab_read objects
		#each Ab_read object contains the V and J segments(s) tha VDJ fasta aligned,
		#the cdr3 sequence and the number of mutations from
		#the key to the dict is the sequence identifier that vdjFasta uses

		#num_Reads is the total number of reads
		#Reads_split_by_V is a dict of dicts  where all the 
		#Ab_read objects are grouped by the V germline 
		#Format:  Reads_split_by_V = {Vgermline :{ID: AB_read object}}


		print "loading sequences..."
		self.Reads = parse_v_j(filepath, num, ABtype)
		self.num_Reads = len(self.Reads)
		print "finding V segments..."
		self.Reads_split_by_V = v_split(self.Reads)

		print "Calculating V segment usage..."
		self.V_freqs = {}
		self.V_fractions = {}
		for Vgerm in self.Reads_split_by_V:
			self.V_freqs[Vgerm] = len(self.Reads_split_by_V[Vgerm])
			self.V_fractions[Vgerm] = len(self.Reads_split_by_V[Vgerm])/float(self.num_Reads)


	def split_V_4(self):
		#used to get counts of V segment usage, just using first number (i.e. IVH4 instead of IGH4-32)
		self.Reads_split_by_f4V = v_first4_split(self.Reads)

		self.f4V_freqs = {}
		self.f4V_fractions = {}
		for Vgerm in self.Reads_split_by_f4V:
			self.f4V_freqs[Vgerm] = len(self.Reads_split_by_f4V[Vgerm])
			self.f4V_fractions[Vgerm] = len(self.Reads_split_by_f4V[Vgerm])/float(self.num_Reads)


	def find_clones(self):
		#This creates a dict of clone objects from all the cdr3's in reads
		#a clone object combines cdr3's that are identical as well as
		#cdr3's within a Hamming distance of 1 from each other (to 
		#account for sequencing errors - see von Budingen, 2012	)
		
		#the keys to the dict are numbers that correspond to the read rank of
		#the clone (i.e. the clone that appears the most in the repertoire has the key '1')
		#the object has the consensus cddr3 sequence, V, J, number of reads and percent of rep
		#that this clone represents

		clone_index = 1
		Vgerm_complete=0
		All_clones = {}
		Clones_split_by_V = {}


		for Vgerm in self.Reads_split_by_V:
			print "finding " + Vgerm + " clones from %s reads" % len(self.Reads_split_by_V[Vgerm])

			

			# Clones_split_by_V[Vgerm] = {}

			# #load in cdr3 sequences from Vgerm
			# all_cdr3s, cdr3_dict = load_cdrs(self.Reads_split_by_V[Vgerm], 'cdr3')
			# all_cdr2s, cdr2_dict = load_cdrs(self.Reads_split_by_V[Vgerm], 'cdr2')
			# all_cdr1s, cdr1_dict = load_cdrs(self.Reads_split_by_V[Vgerm], 'cdr1')
			# all_IDs, ID_dict = load_cdrs(self.Reads_split_by_V[Vgerm], 'ID')

			# #if VDJ didn't find the CDR3, skip to next
			# if all_cdr3s == []:
			# 	Vgerm_complete +=1	
			# 	print str(Vgerm_complete) +'/' +str(len(self.Reads_split_by_V)) + " germlines done!"
			# 	continue

			# #cluster all cdrs by maximum Hamming length = 1 in each
			# T = cluster_into_clones(all_cdr1s, all_cdr2s, all_cdr3s)

			# #list of numbers of clones
			# clones = [x+1 for x in range(np.amax(T))]

			# #find properties of the clone
			# for clone_num in clones: #clones is list of numbers from T

			# 	J, final_seq_1, final_seq_2, final_seq_3, \
			# 	num_reads, percent_reads, ABtype, IDs, Vmut, Jmut, sh = find_clone_props(all_cdr1s,\
			# 															all_cdr2s,\
			# 															all_cdr3s, \
			# 															cdr3_dict,\
			# 															T, \
			# 															self.Reads_split_by_V[Vgerm], \
			# 															self.num_Reads,\
			# 															clone_num,\
			# 															all_IDs)
			# 	#plug it all into a clone object
			# 	All_clones[clone_index] = Clone(V=Vgerm, \
			# 	                               J=J, \
			# 	                               cdr1 = Seq(final_seq_1, generic_protein), \
			# 	                               cdr2 = Seq(final_seq_2, generic_protein), \
			# 	                               cdr3 = Seq(final_seq_3, generic_protein), \
			# 	                               num_reads=num_reads, \
			# 	                               percent_reads=percent_reads, \
			# 	                               ABtype = ABtype,\
			# 	                               IDs = IDs,\
			# 	                               Vmut = Vmut,\
			# 	                               Jmut = Jmut,\
			# 	                               sh = sh)


			# 	Clones_split_by_V[Vgerm][clone_num] = Clone(V=Vgerm, \
			# 	                               			J=J, \
			# 	                               			cdr1 = Seq(final_seq_1, generic_protein), \
			# 	                               			cdr2 = Seq(final_seq_2, generic_protein), \
			# 	                               			cdr3 = Seq(final_seq_3, generic_protein), \
			# 	                               			num_reads=num_reads, \
			# 	                               			percent_reads=percent_reads,\
			# 	                               			ABtype = ABtype,\
			# 	                               			IDs = IDs,\
			# 	                               			Vmut = Vmut,\
			# 	                               			Jmut = Jmut,\
			# 	                               			sh = sh)
			# 	clone_index+=1
			
			# #print progress 
			# Vgerm_complete +=1	
			# print "%s unique clones found " % np.amax(T)
			# print str(Vgerm_complete) +'/' +str(len(self.Reads_split_by_V)) + " germlines done!\n"


		#Order most frequent to least frequent
		All_clones_sorted = order_clones(All_clones)
		Clones_split_by_V_sorted = {}
		for V in Clones_split_by_V:
			Clones_split_by_V_sorted[V] = order_clones(Clones_split_by_V[V])
		
		self.Clones = All_clones_sorted
		self.Clones_split_by_V = Clones_split_by_V_sorted
		self.num_clones = len(All_clones_sorted)

		#find some frequency properties and distributions

		#find the fraction of the repertoir that top clones take up in 1% intervals
		#useful for plotting distributions
		self.clone_distribution = find_distribution(self.Clones, 100, self.num_clones, self.num_Reads)

		#find the fractions of repertoire the top 1, 10 and 100 clones take up 
		self.top_clone_fraction, self.top_10_fraction, self.top_100_fraction = top_clone_fractions(self.Clones, self.num_Reads)

		#find the fraction of the repertoire devoted to each AB type
		#a dict with the keys "IGHG", "IGHM" etc.
		self.ABtype_fractions, self.ABtype_unique = ABtype_fractions(self.Clones, self.num_clones, self.num_Reads)

		clone_split_by_V_fractions = {}




	def find_clusters(self):
		cluster_index = 1
		Vgerm_complete=0
		All_clusters = {}
		read_inds_sum=0
		for Vgerm in self.Clones_split_by_V:
			print "finding " + Vgerm + " clusters from %s clones" % len(self.Clones_split_by_V[Vgerm])

			all_cdr3s, cdr3_dict = load_cdrs(self.Clones_split_by_V[Vgerm], 'cdr3')
			all_cdr2s, cdr2_dict = load_cdrs(self.Clones_split_by_V[Vgerm], 'cdr2')
			all_cdr1s, cdr1_dict = load_cdrs(self.Clones_split_by_V[Vgerm], 'cdr1')
			all_IDs, ID_dict = load_cdrs(self.Clones_split_by_V[Vgerm], 'IDs')

			cdrs=[]
			for i in range(len(all_cdr3s)):
				cdrs.append([all_cdr1s[i], all_cdr2s[i], all_cdr3s[i]])

			T = similarity_cluster(cdrs)



			clusters = [x+1 for x in range(np.amax(T))]

			#find properties of the clone
			for cluster_num in clusters: #clones is list of numbers from T
				read_inds = np.where(T==cluster_num)[0]
				cluster_cdr1s=[]
				cluster_cdr2s=[]
				cluster_cdr3s=[]
				cluster_IDs = []

				for i in read_inds:

					cluster_cdr1s.append(all_cdr1s[i])
					cluster_cdr2s.append(all_cdr2s[i])
					cluster_cdr3s.append(all_cdr3s[i])
					cluster_IDs.append(all_IDs[i])

				Js, ABtypes, num_reads, percent_reads = find_cluster_props(cluster_cdr3s,\
																		cdr3_dict,\
																		self.Clones_split_by_V[Vgerm], \
																		read_inds)


				All_clusters[cluster_index] = Cluster(V=Vgerm, \
				                               Js=Js, \
				                               cdr1s = cluster_cdr1s, \
				                               cdr2s = cluster_cdr2s, \
				                               cdr3s = cluster_cdr3s, \
				                               num_reads=num_reads, \
				                               percent_reads=percent_reads, \
				                               ABtypes = ABtypes,\
				                               IDs = cluster_IDs)

				cluster_index+=1

			
			#print progress 
			Vgerm_complete +=1
			print "%s clusters created " % np.amax(T)	
			print str(Vgerm_complete) +'/' +str(len(self.Reads_split_by_V)) + " germlines done!\n"


		#Order most frequent to least frequent
		All_clusters_sorted = order_clones(All_clusters)		
		self.Clusters = All_clusters_sorted
		self.num_clusters = len(All_clusters_sorted)





	def sh_props(self):
		if not self.Clones:
			print "Run .find_clones() method first"
			return ""

		sh_dict = {'all classes' : [], 'IGHM': [], 'IGHG': [], 'IGHA': [], 'IGHE': [], 'IGHD':[]}

		for clone in self.Clones:
			sh_dict['all classes'].append(self.Clones[clone].sh)
			if self.Clones[clone].ABtype:
				sh_dict[self.Clones[clone].ABtype].append(self.Clones[clone].sh)

		self.sh_dict = sh_dict
		self.median_sh = {'all' : np.median(sh_dict['all classes']),
						'IGHM' : np.median(sh_dict['IGHM']), 
						'IGHG' : np.median(sh_dict['IGHG']),
						'IGHA' : np.median(sh_dict['IGHA']),
						'IGHE' : np.median(sh_dict['IGHE']),
						'IGHD' : np.median(sh_dict['IGHD']),}


		def sh_plots(data, label, color):
			if not data:
				print "no %s found"%label
				return ""
			bins = np.linspace(0, max(data), max(data)+1)

			plt.rc('xtick', labelsize=15) 
			plt.rc('ytick', labelsize=15)
			# plt.rc('legend',fontsize=15)
			# plt.rc('axes',labelsize=15)

			plt.figure(num=None, figsize=(7, 5), dpi=80)
			plt.hist(np.array(data), bins, color=color, alpha=.5, histtype = 'stepfilled')

			plt.xlabel("Somatic Hypermutations", fontsize=15)
			plt.ylabel("number of clones", fontsize=15) 
			plt.legend([label],  fontsize=15)
			plt.title(label, fontsize = 20)

			plt.show()


		#Make histograms for all classes
		try:
			c1, c2, c3, c4, c5, c6 = sns.color_palette("Set1", 6)
			colors = [c1, c2, c3, c4, c5, c6]
		except:
			colors = ['r','b','g','y','k', 'c']
		
		color_index = 0
		for key in list(sh_dict.keys()):
			sh_plots(sh_dict[key], key, colors[color_index])
			color_index+=1


		#Make box and whisker plot
		box_data = []
		box_labels = []
		for ab_class in ['IGHD', 'IGHM', 'IGHG', 'IGHA','IGHE']: #some problem with IGHE - fix later
			if len(sh_dict[ab_class])>5:
				box_data.append(np.array(sh_dict[ab_class]))
				box_labels.append(ab_class)


		plt.figure(num=None, figsize=(7, 5), dpi=80)
		plt.rc('xtick', labelsize=15) 
		plt.rc('ytick', labelsize=15)

		plt.boxplot(box_data)

		plt.xticks([x+1 for x in range(len(box_labels))], box_labels,  fontsize=15)
		plt.ylabel('Somatic Hypermutations', fontsize=15)



	def plot_V_freqs(self, f4=False):
		if f4:
			data = self.f4V_freqs
		else:
			data = self.V_freqs

		plt.bar(range(len(data)), data.values())
		plt.xticks(range(len(data)), data.keys(), rotation=90)
		plt.xlabel('V-Germlines')
		plt.ylabel('Reads')
		plt.show()




	def plot_V_fractions(self, f4=False):
		if f4:
			data = self.f4V_fractions
		else:
			data = self.V_fractions

		plt.bar(range(len(data)), data.values())
		plt.xticks(range(len(data)), data.keys(), rotation=90)
		plt.xlabel('V-Germlines')
		plt.ylabel('Fraction of Repertoire')
		plt.show()		








