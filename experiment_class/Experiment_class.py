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

		#Creates a color dictionary that fades per sample and is based on the input color
		self.color_dict = {}
		for group in self.groups:
		    color = sns.color_palette(group['color'], len(group['samples']))
		    for i, sample in enumerate(group['samples']):
		        self.color_dict[sample] = color[i]

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
		        

	def violins(self): #SHIFT THIS TO DATAFRAME and change the coloring on this

		#Creates a dict of dicts with the following structure SHM_dict[antibody class][group name][sample name] and the resulting 
		#object is a list with the number of somatic hyper mutations 
		antibody_class_list = ['all classes', 'IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD']

		SHM_dict = {}
		for key in antibody_class_list:  
		    SHM_dict[key] = {}
		    for igroup,group in enumerate(self.groups):
		        SHM_dict[key][group['name']] = {} 
		        for sample in group['samples']:
		            SHM_dict[key][group['name']][sample] = self.groups[igroup]['sample data'][sample]['sh_dict'][key]

		self.SHM_dict = SHM_dict

		bins = np.linspace(0, 50, 50+1)

		plot_input = {}
		for key in self.SHM_dict: 
			plot_input[key] = []
			no_clones = []
			x_labels = []
			colors = []
			for group in self.SHM_dict[key]:
				for sample in self.SHM_dict[key][group]:
		            
					data = self.SHM_dict[key][group][sample]

					if len(data)<4: 
						no_clones.append(sample)
					else:
						plot_input[key].append(data)
						x_labels.append(sample)
						colors.append(self.color_dict[sample])
			
			#print key, plot_input[key]
			if plot_input[key]:
				plt.title(key)
				sns.violinplot(plot_input[key], color=colors, names=x_labels)
				plt.show()
			if no_clones:
				print "There wasn't enough data for %s in %s" % (no_clones, key)

	#def SHM_density(self):

	def clone_bars(self):

		self.clone_DF.plot(x='Sample', kind='bar', figsize=(20,6));

	def clone_cum_dist(self):

		colors = sns.color_palette("hls", len(self.groups))
		patches = []
		for igroup, group in enumerate(self.groups):
		    patches.append(mpatches.Patch(color=colors[igroup], label=group['name']))
		    for sample in group['samples']:
		        data = self.groups[igroup]['sample data'][sample]['clone_distribution']
		        plt.plot(np.array(data).cumsum(), color=colors[igroup])
		plt.legend(handles=patches, loc=4);


	def split_SHM_DF(self):
		SHM_columns = ['Antibody', 'Group', 'Sample', 'Data']
		antibody_class_list = ['all classes', 'IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD']

		ls_split = []
		for key in antibody_class_list:
			for igroup,group in enumerate(self.groups):
				for sample in group['samples']:
					for i in self.groups[igroup]['sample data'][sample]['sh_dict'][key]:
						ls_split.append([key, group['name'], sample, i])
		DF_split = pd.DataFrame(ls_split, columns = SHM_columns)

		self.SHM_DF_split = DF_split

	def SHM_factor_plots(self):

		try: sns.factorplot('Sample', y='Data', hue='Group', data=self.SHM_DF_split, col='Antibody'); 	
		except: print "run self.split_SHM_DF"

	def pruned_vs_full(self):

		colors = sns.color_palette("hls", len(self.groups))
		plt.figure(figsize=(20,6))

		for igroup, group in enumerate(self.groups):
		    for sample in group['samples']:
		        fullDF = pd.DataFrame.from_dict(self.groups[igroup]['sample data'][sample]['Full_Tree_Size'])
		        prunedDF = pd.DataFrame.from_dict(self.groups[igroup]['sample data'][sample]['Pruned_Tree_Size'])
		        
		        prune = prunedDF.values.reshape(702)
		        full = fullDF.values.reshape(702)
		        
		        plt.scatter(prune, full, s=10, c=colors[igroup]);
		        
		plt.show()

	def pruned_vs_total_reads(self):

		colors = sns.color_palette("hls", len(self.groups))
		plt.figure(figsize=(20,6))

		for igroup, group in enumerate(self.groups):
		    for sample in group['samples']:
		        fullDF = pd.DataFrame.from_dict(self.groups[igroup]['sample data'][sample]['Full_Tree_Size'])
		        #prunedDF = pd.DataFrame.from_dict(test.groups[igroup]['sample data'][sample]['Pruned_Tree_Size'])
		        
		        tree_size_list = []
		        total_reads = []
		        for key in self.groups[igroup]['sample data'][sample]['VJ_freqs']:
		            V = key.split('_')[0]
		            J = key.split('_')[1]

		            tree_size_list.append(fullDF[J].loc[V])
		            total_reads.append(test.groups[igroup]['sample data'][sample]['VJ_freqs'][key])
		            
		        plt.scatter(tree_size_list, total_reads, s=10, c=colors[igroup]);
		        
		plt.show()
		        
        
