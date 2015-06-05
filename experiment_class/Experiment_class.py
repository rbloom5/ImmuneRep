import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd 
from boto.s3.connection import S3Connection
from boto.s3.key import Key

class Experiment:
	def __init__(self,groups):

		#Opens up an S3 connection and sets buckets
		#dont share the strings below with anyone!
		conn = S3Connection('AKIAJ2TEUHQV2LHU7XQQ','VKzoYINBlvZi5uiAIeAnJG5fgLedQPFrmMCpSfBp')
		json_buck = conn.get_bucket('rep-seq-jsons')
		info_buck = conn.get_bucket('info-docs')

		#Checks for groups with no input samples and looks for a default
		info_buck.get_key("ExperimentClassDefaultGroups.txt").get_contents_to_filename("ExperimentClassDefaultGroups.txt")
		with open("ExperimentClassDefaultGroups.txt") as default_groups_location:
			defaultgroups=json.load(default_groups_location) 

		for igroup, group in enumerate(groups):
			if 'samples' not in group:
				groups[igroup]['samples']=defaultgroups[group['name']]

		#Downloads all the JSONs from AWS to the current directory
		print "downloading jsons"
		for key in json_buck.list():
			key.get_contents_to_filename('tempJSONs/'+key.name)

		#Pulls in the sample data for each sample
		for igroup, group in enumerate(groups):
			groups[igroup]['sample data'] = {}
			for sample in group['samples']:
				groups[igroup]['sample data'][sample]=json.load(open('tempJSONs/'+sample+'.json'))

		self.groups = groups

		#Creates a Pandas DataFrame with the columns [Antibody  Group  Sample  Data], where the index is the sample name and Data is a list of
		#number of somatic hypermutations
		SHM_columns = ['Antibody', 'Group', 'Sample', 'Data']
		antibody_class_list = ['all classes', 'IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD']

		ls = []
		for key in antibody_class_list:
			for igroup,group in enumerate(self.groups):
				for sample in group['samples']:
					SHM_data = self.groups[igroup]['sample data'][sample]['sh_dict'][key]
					ls.append([key, group['name'], sample, SHM_data])

		self.SHM_DF = pd.DataFrame(ls, columns = SHM_columns)

		#self.SHM_DF.set_index('Sample')

		#Creates a Panda DataFrame for the clone information. The columns are [Group  Sample  Top 1  Top 10  Top 100]
		Clone_columns = ['Group', 'Sample', 'Top 1', 'Top 10', 'Top 100']

		clone_ls = []
		for igroup, group in enumerate(self.groups):
		    for sample in group['samples']:
		        top = self.groups[igroup]['sample data'][sample]['top_clone_fraction']
		        top_10 = self.groups[igroup]['sample data'][sample]['top_10_fraction']
		        top_100 = self.groups[igroup]['sample data'][sample]['top_100_fraction']
		        
		        clone_ls.append([group['name'], sample, top, top_10, top_100])

		self.clone_DF = pd.DataFrame(clone_ls, columns=Clone_columns)

#Investigating clone diversity		

	def clone_bars(self):

		self.clone_DF.plot(x='Sample', kind='bar', figsize=(20,6));

	def clone_cum_dist(self):

		plt.figure(figsize=(20,6))

		colors = sns.color_palette("hls", len(self.groups))
		patches = []
		for igroup, group in enumerate(self.groups):
		    patches.append(mpatches.Patch(color=colors[igroup], label=group['name']))
		    for sample in group['samples']:
		        data = self.groups[igroup]['sample data'][sample]['clone_distribution']
		        plt.plot(np.array(data).cumsum(), color=colors[igroup])
		plt.legend(handles=patches, loc=4);

#Investigating Somatic Hyper Mutation	

	def class_SHM(self):

		SHM_columns = ['Antibody', 'Group', 'Sample', 'Data']
		antibody_class_list = ['all classes', 'IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD']

		ls_split = []
		for key in antibody_class_list:
			for igroup,group in enumerate(self.groups):
				for sample in group['samples']:
					for i in self.groups[igroup]['sample data'][sample]['sh_dict'][key]:
						ls_split.append([key, group['name'], sample, i])
		self.SHM_DF_split = pd.DataFrame(ls_split, columns = SHM_columns)

		sns.factorplot('Group', 'Data', data=self.SHM_DF_split, hue='Antibody' , kind='box', size = 8, aspect=2)

	def class_fractions(self, antibody_type=False):
		columns = ['Group', 'Sample', 'Antibody', 'Data']
		antibody_class_list = ['IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD']

		ls = []
		for key in antibody_class_list:
		    for igroup,group in enumerate(self.groups):
		        for sample in group['samples']:
		            data = self.groups[igroup]['sample data'][sample]['ABtype_fractions'][key]
		            ls.append([group['name'], sample, key, data])
		            
		self.class_fraction_DF = pd.DataFrame(ls, columns = columns)

		if antibody_type:
			isolate_class = self.class_fraction_DF[self.class_fraction_DF.Antibody == antibody_type]
			sns.factorplot('Sample', 'Data', data=isolate_class, hue='Group', kind='bar', size = 8, aspect=2);
		else:
			sns.factorplot('Group', 'Data', data=self.class_fraction_DF, hue='Antibody', kind='bar', size = 8, aspect=2);

		#df = self.class_fraction_DF[self.class_fraction_DF.Antibody == 'IGHM']

#Tree analysis

	def tree_scatter(self, x, y, size_variable=False):
		#Sets the color palette from Seaborn
	    colors = sns.color_palette("hls", len(self.groups))
	    plt.figure(figsize=(20,10))

	    columns = ['VJ', 'VJ Freq', 'Tree Size', 'Pruned Tree Size', 'VJ Fraction', 'Generations', 'Avg Mutations', 'STD Mutations']
	    patches = []

	    for igroup, group in enumerate(self.groups):
	    	#assigns each group it's own color from the color pallete
	        patches.append(mpatches.Patch(color=colors[igroup], label=group['name']))

	        for sample in group['samples']:

	        	#Turns the tree dicts into DataFrames
	            fullDF = pd.DataFrame.from_dict(self.groups[igroup]['sample data'][sample]['Full_Tree_Size'])
	            prunedDF = pd.DataFrame.from_dict(self.groups[igroup]['sample data'][sample]['Pruned_Tree_Size'])

	            #Iterates through every VJ pair in the VJ generations dict. VJ generations has the smallest subset of pairs becuase it is an 
	            #output of Immunitree which cannot run on VJ pairs with only one clone. Therefore the pairs need to be filtered down. 
	            data = []
	            for key in self.groups[igroup]['sample data'][sample]['VJ_generations']:
	                V = key.split('_')[0]
	                J = key.split('_')[1]

	                VJ_freq = self.groups[igroup]['sample data'][sample]['VJ_freqs'][key]
	                VJ_frac = self.groups[igroup]['sample data'][sample]['VJ_fractions'][key]
	                tree_size = fullDF[J].loc[V]
	                pruned_tree_size = prunedDF[J].loc[V]
	                generations = self.groups[igroup]['sample data'][sample]['VJ_generations'][key]
	                avg_mutations = self.groups[igroup]['sample data'][sample]['VJ_shm_mean'][key]
	                std_mutations = self.groups[igroup]['sample data'][sample]['VJ_shm_stdev'][key]
	                    
	                data.append([key, VJ_freq, tree_size, pruned_tree_size, VJ_frac*10000, generations, avg_mutations, std_mutations])
	                    
	            plotting_DF = pd.DataFrame(data, columns = columns).set_index('VJ')

	            #Looks to see if user wants circle size to repersent something
	            if size_variable != False : circle_size = plotting_DF[size_variable].tolist()
	            else: circle_size = 100

	            plt.scatter(plotting_DF[x].tolist(), plotting_DF[y].tolist(), s=circle_size, c=colors[igroup], alpha=0.5);

	    #makes it so that 0,0 is the lower left corner
	    plt.xlim(0)
	    plt.ylim(0)
	    plt.ylabel(y)
	    plt.xlabel(x)

	    #Creates a legend with patches
	    plt.legend(handles=patches, loc=4)
	    plt.show()
			        
	        
