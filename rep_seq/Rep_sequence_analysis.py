#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import re
import os
import numpy as np
from collections import Counter
import itertools
from joblib import Parallel, delayed  
import multiprocessing
from multiprocessing import Pool, freeze_support
from collections import OrderedDict
from ete2 import Tree
import glob
import operator



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
# import ab_classes 
# reload(ab_classes)
from ab_classes import *
      

class Rep_seq:

	def __init__(self, filepath, num = False, ABtype = []): #filepath must be file with extension XXX.VDJ.H3.L3.CH1.fa from VDJfasta
		#parses the data from XXX.VDJ.H3.L3.CH1.fa into usable format
		#num is an optional input where you can set the max number of reads to take
		#from the file if you don't want them all

		#self.Reads is a dict of Ab_read objects
		#each Ab_read object contains the V and J segments(s) tha VDJ fasta aligned,
		#the cdr3 sequence and the number of mutations 
		#the key to the dict is the sequence identifier that vdjFasta uses

		#num_Reads is the total number of reads
		#Reads_split_by_VJ is a dict of dicts  where all the 
		#Ab_read objects are grouped by the VJ germlines 
		#Format:  Reads_split_by_VJ = {VJgermline :{ID: AB_read object}}

		self.features={}

		self.filepath = filepath
		print "loading sequences..."
		self.Reads = parse_v_j(filepath, num, ABtype)
		self.features['num_Reads'] = len(self.Reads)
		print "finding V segments..."
		self.Reads_split_by_VJ = vj_split(self.Reads)

		print "Calculating V segment usage..."
		self.features['VJ_freqs'] = {}
		self.features['VJ_fractions'] = {}
		for germ in self.Reads_split_by_VJ:
			self.features['VJ_freqs'][germ] = len(self.Reads_split_by_VJ[germ])
			self.features['VJ_fractions'][germ] = len(self.Reads_split_by_VJ[germ])/float(self.features['num_Reads'])


	def split_V_4(self):
		#used to get counts of V segment usage, just using first number (i.e. IVH4 instead of IGH4-32)
		self.Reads_split_by_f4V = v_first4_split(self.Reads)

		self.f4V_freqs = {}
		self.features['f4VJ_fractions'] = {}
		for germ in self.Reads_split_by_f4V:
			self.f4V_freqs[germ] = len(self.Reads_split_by_f4V[germ])
			self.features['f4VJ_fractions'][germ] = len(self.Reads_split_by_f4V[germ])/float(self.features['num_Reads'])


	def find_clones(self,parallel=False):
		#This creates a dict of clone objects from all the sequences in reads
		#a clone object combines sequences that are identical as well as
		#within a Hamming distance of 1 from each other (to 
		#account for sequencing errors - see von Budingen, 2012	)
		
		#the keys to the dict are numbers that correspond to the read rank of
		#the clone (i.e. the clone that appears the most in the repertoire has the key '1')
		#the object has the consensus cddr3 sequence, V, J, number of reads and percent of rep
		#that this clone represents
		print "Finding clones..."
		clone_index = 1
		germ_complete=0
		All_clones = {}
		# All_clones = []
		Clones_split_by_VJ = {}


		if parallel:
			num_cores = multiprocessing.cpu_count()/4
			pool = Pool(processes=num_cores)


			results = pool.map(reads_to_clones, itertools.izip(itertools.repeat(self.Reads_split_by_VJ),\
															list(self.Reads_split_by_VJ.keys()), \
															itertools.repeat(self.features['num_Reads'])))
			All_clones = [result[0] for result in results]
			ind_Clones_split_by_VJ = [result[1] for result in results]		
			All_clones = list(itertools.chain.from_iterable(All_clones))
			All_clones = {key+1: value for key, value in enumerate(All_clones)}

			for item in ind_Clones_split_by_VJ:
				Clones_split_by_VJ[item.keys()[0]] = item[item.keys()[0]]

			pool.close()
			pool.join()


		else:
			for germ in self.Reads_split_by_VJ:
				print "finding " + germ + " clones from %s reads" % len(self.Reads_split_by_VJ[germ])


				Clones_split_by_VJ[germ] = {}

				#load in cdr3 sequences from germ
				all_cdr3s, cdr3_dict = load_cdrs(self.Reads_split_by_VJ[germ], 'cdr3')
				all_cdr2s, cdr2_dict = load_cdrs(self.Reads_split_by_VJ[germ], 'cdr2')
				all_cdr1s, cdr1_dict = load_cdrs(self.Reads_split_by_VJ[germ], 'cdr1')
				all_IDs, ID_dict = load_cdrs(self.Reads_split_by_VJ[germ], 'ID')

				#if VDJ didn't find the CDR3, skip to next
				if all_cdr3s == []:
					germ_complete +=1	
					print str(germ_complete) +'/' +str(len(self.Reads_split_by_VJ)) + " germlines done!"
					continue

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
																			self.Reads_split_by_VJ[germ], \
																			self.features['num_Reads'],\
																			clone_num,\
																			all_IDs)
					#plug it all into a clone object
					All_clones[clone_index] = Clone(V=germ, \
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
					                               sh = sh )


					Clones_split_by_VJ[germ][clone_num] = Clone(V=germ, \
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
					                               			sh = sh )
					clone_index+=1
				
				#print progress 
				germ_complete +=1	
				print "%s unique clones found " % np.amax(T)
				print str(germ_complete) +'/' +str(len(self.Reads_split_by_VJ)) + " germlines done!\n"


		#Order most frequent to least frequent
		All_clones_sorted = order_clones(All_clones)
		Clones_split_by_VJ_sorted = {}
		for V in Clones_split_by_VJ:
			Clones_split_by_VJ_sorted[V] = order_clones(Clones_split_by_VJ[V])
		
		self.Clones = All_clones_sorted
		self.Clones_split_by_VJ = Clones_split_by_VJ_sorted
		self.features['num_clones'] = len(All_clones_sorted)




	def tree(self):
		#a function that will run immunitree to make lineage trees
		#the process is to first make a file for each VJ pair then write the
		# sequences for all the clones from that VJ pair to the file
		# then shoot all the files into a compiled version of immunitree
		# get the output, then delete all the files


		#make somewhat descriptive directory and file names
		head, tail = os.path.split(self.filepath[0])
		dirstring = tail.split('.')[0]+'_vj_files'
		os.system("mkdir "+dirstring)
		os.system('mkdir /home/ubuntu/tree_output')
		os.system('mkdir /home/ubuntu/parsed_fasta')

		file_list = []
		# make fasta file for all clones in each VJ pair
		print "writing files for tree creation"
		for germ in self.Clones_split_by_VJ:
			filestring = dirstring+'/immTree_'+germ+'.fasta'
			file_list.append(filestring)
			with open(filestring, 'w') as f:
				for clone in self.Clones_split_by_VJ[germ]:
					best_id = find_best_id(self.Clones_split_by_VJ[germ][clone], self.Reads_split_by_VJ[germ])
					find_and_write(best_id, f, self.filepath[0])


		# set up a parallell processing pool
		num_cores = multiprocessing.cpu_count()
		pool = Pool(processes=num_cores)
		nIter = 300 #should do 300+ for large repertoires
		fstrings_for_pool=[]

		matlab_call = '/home/ubuntu/imm_tree_test/run_run_immunitree.sh \
						/usr/local/v83 '
						#/home/ubuntu/SRR1383448_vj_files/immTree_IGHV1-2_IGHJ4.fasta 50'

		#make list of all the bash calls (one for each VJ file we input to immunitree)
		for f in file_list:
			fstrings_for_pool.append(matlab_call + f + ' %s'%nIter)

		#run immunitree on all cores, the .node files are stored in S3 
		print "Running Immunitree.  This may take a while..."
		try:
			results = pool.map(os.system, fstrings_for_pool)
		except:
			pass

		pool.close()
		pool.join()
		#delete all the temp vj-files
		os.system('sudo rm -f -r %s'%dirstring)

		# This section creates a dict where the key is the VJ pair (seperated by _) and the object has the 
		# ete2 tree structure. Each node in the tree structure has associated with it; a name, size, mutations from parent, and
		# a SEQ object. The second dict is identical to the first, it just has the leaves removed
		self.tree_dict = {}
		self.pruned_tree_dict ={}

		for entry in glob.glob('/home/ubuntu/tree_output/*.node.txt'):
			try:
				VJ_name = entry.split('/')[-1].split('.')[0][8:]
				self.tree_dict[VJ_name] = convert_immunitree_to_ete2(entry)
				self.pruned_tree_dict[VJ_name] = convert_immunitree_to_ete2(entry)
				for leaf in self.pruned_tree_dict[VJ_name].get_leaves(): leaf.detach()
			except:
				continue



	def find_clusters(self):
		cluster_index = 1
		germ_complete=0
		All_clusters = {}
		read_inds_sum=0
		for germ in self.Clones_split_by_VJ:
			print "finding " + germ + " clusters from %s clones" % len(self.Clones_split_by_VJ[germ])

			all_cdr3s, cdr3_dict = load_cdrs(self.Clones_split_by_VJ[germ], 'cdr3')
			all_cdr2s, cdr2_dict = load_cdrs(self.Clones_split_by_VJ[germ], 'cdr2')
			all_cdr1s, cdr1_dict = load_cdrs(self.Clones_split_by_VJ[germ], 'cdr1')
			all_IDs, ID_dict = load_cdrs(self.Clones_split_by_VJ[germ], 'IDs')

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
																		self.Clones_split_by_VJ[germ], \
																		read_inds)


				All_clusters[cluster_index] = Cluster(V=germ, \
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
			germ_complete +=1
			print "%s clusters created " % np.amax(T)	
			print str(germ_complete) +'/' +str(len(self.Reads_split_by_VJ)) + " germlines done!\n"


		#Order most frequent to least frequent
		All_clusters_sorted = order_clones(All_clusters)		
		self.Clusters = All_clusters_sorted
		self.features['num_clusters'] = len(All_clusters_sorted)





	def get_stats(self):
		if not self.Clones:
			print "Run .find_clones() method first"
			return ""

		#### Get Somatic Hypermutation Stats ######
		sh_dict = {'all classes' : [], 'IGHM': [], 'IGHG': [], 'IGHA': [], 'IGHE': [], 'IGHD':[]}

		for clone in self.Clones:  #is reads or clones better?
			sh_dict['all classes'].append(self.Clones[clone].sh)
			if self.Clones[clone].ABtype and self.Clones[clone].ABtype in sh_dict:
				sh_dict[self.Clones[clone].ABtype].append(self.Clones[clone].sh)

		self.features['sh_dict'] = sh_dict
		self.features['median_sh'] = {'all' : np.median(sh_dict['all classes']),
						'IGHM' : np.median(sh_dict['IGHM']), 
						'IGHG' : np.median(sh_dict['IGHG']),
						'IGHA' : np.median(sh_dict['IGHA']),
						'IGHE' : np.median(sh_dict['IGHE']),
						'IGHD' : np.median(sh_dict['IGHD']),}

		

		#### Frequency properties and distributions #####

		#find the fraction of the repertoir that top clones take up in 1% intervals
		#useful for plotting distributions
		self.features['clone_distribution'] = find_distribution(self.Clones, 100, self.features['num_clones'], self.features['num_Reads'])

		#find the fractions of repertoire the top 1, 10 and 100 clones take up 
		self.features['top_clone_fraction'], self.features['top_10_fraction'], self.features['top_100_fraction'] = top_clone_fractions(self.Clones, self.features['num_Reads'])

		#find the fraction of the repertoire devoted to each AB type
		#a dict with the keys "IGHG", "IGHM" etc.
		self.features['ABtype_fractions'], self.features['ABtype_unique'] = ABtype_fractions(self.Clones, self.features['num_clones'], self.features['num_Reads'])




	def plots(self):

		plt.rc('xtick', labelsize=15) 
		plt.rc('ytick', labelsize=15)



		#### Make histograms for all classes #####
		try:
			c1, c2, c3, c4, c5, c6 = sns.color_palette("Set1", 6)
			colors = [c1, c2, c3, c4, c5, c6]
		except:
			colors = ['r','b','g','y','k', 'c']
		
		color_index = 0
		for key in list(self.features['sh_dict'].keys()):
			sh_plots(self.features['sh_dict'][key], key, colors[color_index])
			color_index+=1




		###### Make box and whisker plot ########
		box_data = []
		box_labels = []
		for ab_class in ['IGHD', 'IGHM', 'IGHG', 'IGHA','IGHE']: 
			if len(self.features['sh_dict'][ab_class])>5:
				box_data.append(np.array(self.features['sh_dict'][ab_class]))
				box_labels.append(ab_class)

		plt.figure(num=None, figsize=(7, 5), dpi=80)
		plt.rc('xtick', labelsize=15) 
		plt.rc('ytick', labelsize=15)
		plt.boxplot(box_data)
		plt.xticks([x+1 for x in range(len(box_labels))], box_labels,  fontsize=15)
		plt.ylabel('Somatic Hypermutations', fontsize=15)




		######## Plot Ab type fractions #############
		bar_plt(self.features['ABtype_fractions'], 'Fraction of reads', 'AB type Read Fractions')
		bar_plt(self.features['ABtype_unique'], 'Fraction of unique clones', 'AB type Fraction of Unique Clones', color="#F08080")




		######## Plot clone fraction distributions ##########

		clone_fractions = OrderedDict([('top clone', self.features['top_clone_fraction']), \
							('top 10 clones', self.features['top_10_fraction']), \
							('top 100 clones', self.features['top_100_fraction'])])

		bar_plt(clone_fractions, 'Fraction of total reads', 'Fractions of Repertoire for Top Clones')



		######## Plot cumulative clone fractions ########
		plt.figure(num=None, figsize=(7, 5), dpi=80)
		plt.plot(np.array(self.features['clone_distribution']).cumsum())
		plt.xlabel('Top Clones (in percentage)', fontsize=15)
		plt.ylabel('Cumulative Fraction of Reads',fontsize=15)



		#######  Plot VJ fractions for 10 most common VJ pairs ##########
		self.plot_VJ_fractions(top_ten=True)
		







	def plot_VJ_freqs(self, f4=False, top_ten=False):
		plt.figure(num=None, figsize=(7, 5), dpi=80)
		plt.rc('xtick', labelsize=15) 
		plt.rc('ytick', labelsize=15)

		if top_ten:
			data = dict(sorted(self.features['VJ_freqs'].iteritems(), key=operator.itemgetter(1), reverse=True)[:10])
		elif f4:
			data = self.features['f4VJ_freqs']
		else:
			data = self.features['VJ_freqs']

		plt.bar(range(len(data)), data.values(), align='center')
		plt.xticks(range(len(data)), data.keys(), rotation=90)
		plt.xlabel('VJ-Germlines', fontsize=15)
		plt.ylabel('Reads',fontsize=15)
		plt.show()






	def plot_VJ_fractions(self, f4=False, top_ten=False):
		plt.figure(num=None, figsize=(7, 5), dpi=80)
		plt.rc('xtick', labelsize=15) 
		plt.rc('ytick', labelsize=15)

		if top_ten:
			data = dict(sorted(self.features['VJ_fractions'].iteritems(), key=operator.itemgetter(1), reverse=True)[:10])
		elif f4:
			data = self.features['f4VJ_fractions']
		else:
			data = self.features['VJ_fractions']

		plt.bar(range(len(data)), data.values(), align='center')
		plt.xticks(range(len(data)), data.keys(), rotation=90)
		plt.xlabel('VJ-Germlines', fontsize=15)
		plt.ylabel('Fraction of Repertoire', fontsize=15)
		plt.show()




	def output_features(self,outpath):
		import json
		outfile=open(outpath,'wb')
		json.dump(self.features,outfile,indent=1,sort_keys=True)









