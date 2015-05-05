#!/usr/bin/python
import json
import pickle

from rep_stat_functions import *

def update_rep_stats(reps = default):
	if reps == default:
		reps = [

		###  I'm thinking this is a list of all the id's of all of our reps

		]


	#run through each rep, load it in, update features dict with any new features, and dump it into a new json

	for rep in reps:
		features_dict = json.load()!!!!!!
		with open(rep, 'wb') as rep:
			Rep = pickle.load(rep)


		## now loop through all of our functions and see if it is in the json
		## if not, then calculate it

		if "Full_Tree_Size" not in features_dict:
			features_dict["Full_Tree_Size"] = calculate_tree_size(Rep, pruned=False) #function stored in rep_stat_functions
		

		if other_feature_name not in features_dict:
			features_dict[other_feature_name] = function_that_calculates_other_feature(Rep)


			#  etc. etc. for all interesting funcitons


		# once we run through everything, dump the updated features dict into a new json
		
		with open(outpath,'wb') as outfile:
			json.dump(self.features_dict,outfile,indent=1,sort_keys=True)