#!/usr/bin/python
import json
import pickle
from boto.s3.connection import S3Connection
from boto.s3.key import Key
import os
import subprocess
import sys
sys.path.append(os.path.dirname(os.getcwd()))
from rep_stats_functions import *


def update_rep_stats(reps = 'default'):

	if reps == 'default':

		#dont share the strings below with anyone!
		print "connecting to S3"
		conn = S3Connection('AKIAJ2TEUHQV2LHU7XQQ','VKzoYINBlvZi5uiAIeAnJG5fgLedQPFrmMCpSfBp')
		json_buck = conn.get_bucket('rep-seq-jsons')
		json_buck_list = conn.get_bucket('rep-seq-jsons').list()
		obj_buck_list = conn.get_bucket('rep-seq-objects').list()
		json_names = [str(key.name) for key in json_buck_list]
	else:
		#dont share the strings below with anyone!
		print "connecting to S3"
		conn = S3Connection('AKIAJ2TEUHQV2LHU7XQQ','VKzoYINBlvZi5uiAIeAnJG5fgLedQPFrmMCpSfBp')
		json_buck = conn.get_bucket('rep-seq-jsons')
		json_buck_list = [key for key in conn.get_bucket('rep-seq-jsons').list() if str(key.name)[:-5] in reps] 
		obj_buck_list = [key for key in conn.get_bucket('rep-seq-objects').list() if str(key.name)[:-4] in reps]
		json_names = [str(key.name) for key in json_buck_list]

	# The jsons are small, so just download them all at once
	print "downloading jsons"
	for key in json_buck_list:
		key.get_contents_to_filename(key.name)


	#run through each rep, download it, update features dict with any new features, and dump it into a new json
	index=0
	for key in obj_buck_list:

		#get rep from s3
		key.get_contents_to_filename(key.name)
		rep = str(key.name)
		print "updating ", rep[:-4]

		#load in the correct features dict and rep
		features_dict = json.load(open(json_names[index]))
		with open(rep) as r:
			Rep = pickle.load(r)

		## now loop through all of our functions and see if it is in the json
		## if not, then calculate it

		if "Full_Tree_Size" not in features_dict:
			features_dict["Full_Tree_Size"] = calculate_tree_size(Rep, pruned=False)


		if "Pruned_Tree_Size" not in features_dict:
			features_dict["Pruned_Tree_Size"] = calculate_tree_size(Rep, pruned=True)

		if "clone_cdr3_lengths" not in features_dict:
			features_dict["clone_cdr3_lengths"] = calculate_cdr_lengths(Rep.Clones)
			features_dict["read_cdr3_lengths"] = calculate_cdr_lengths(Rep.Reads)

		if "VJ_shm_mean" not in features_dict:
			features_dict['VJ_shm_mean'], features_dict['VJ_shm_stdev'] = calculate_vj_shm(Rep.Clones_split_by_VJ)

		if "avg_node_mutation_length" not in features_dict:
			features_dict["avg_node_mutation_length"] = avg_node_mutation_length(Rep.tree_dict)
			
		if "VJ_generations" not in features_dict:
			features_dict["VJ_generations"] = generations(Rep.tree_dict)

		if "leafiness_matrix" not in features_dict:
			features_dict["leafiness_matrix"] = leafiness_matrix(Rep.tree_dict)

		if "leafiness_global" not in features_dict:
			features_dict["leafiness_global"] = leafiness_global(Rep.tree_dict)

		if "generations_matrix" not in features_dict:
			features_dict["generations_matrix"] = generations_matrix(features_dict["VJ_generations"])

		if "generations_global" not in features_dict:
			features_dict["generations_global"] = generations_global(features_dict["VJ_generations"])

		
		

		# if other_feature_name not in features_dict:
		# 	features_dict[other_feature_name] = function_that_calculates_other_feature(Rep)


			#  etc. etc. for all interesting funcitons




		# once we run through everything, dump the updated features dict into a new json
		with open(json_names[index],'w') as outfile:
			json.dump(features_dict,outfile,indent=1,sort_keys=True)
		
		j=Key(json_buck)
		j.name = json_names[index]
		j.set_contents_from_filename(json_names[index])

		os.system('rm %s'%json_names[index])
		os.system('rm %s'%rep)
		index+=1
		print "finished updating %s \n"%rep[:-4]
	conn.close()


