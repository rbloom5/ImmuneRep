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


def fix_trees(reps = 'default'):

	if reps == 'default':

		#dont share the strings below with anyone!
		print "connecting to S3"
		conn = S3Connection('AKIAJ2TEUHQV2LHU7XQQ','VKzoYINBlvZi5uiAIeAnJG5fgLedQPFrmMCpSfBp')
		json_buck = conn.get_bucket('rep-seq-jsons')
		obj_buck = conn.get_bucket('rep-seq-objects')
		json_names = [str(key.name) for key in json_buck.list()]

	# The jsons are small, so just download them all at once
	print "downloading jsons"
	for key in json_buck.list():
		key.get_contents_to_filename(key.name)


	#run through each rep, download it, update features dict with any new features, and dump it into a new json
	index=0
	for key in obj_buck.list():

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

		
		features_dict["Full_Tree_Size"] = calculate_tree_size(Rep, pruned=False)

		features_dict["Pruned_Tree_Size"] = calculate_tree_size(Rep, pruned=True)

	
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


