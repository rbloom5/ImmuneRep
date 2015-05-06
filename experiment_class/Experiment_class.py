import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 
from boto.s3.connection import S3Connection
from boto.s3.key import Key

class Experiment:
	def __init__(self,groups,defaultgroupspath='defaultgroups.txt'):
		defaultgroups=json.load(open(defaultgroupspath)) 

		#Checks for groups with no input samples and looks for a default
		for igroup, group in enumerate(groups):
			if 'samples' not in group:
				groups[igroup]['samples']=defaultgroups[group['name']]

		#Downloads all the JSONs from AWS to the current directory
		#dont share the strings below with anyone!
		conn = S3Connection('AKIAJ2TEUHQV2LHU7XQQ','VKzoYINBlvZi5uiAIeAnJG5fgLedQPFrmMCpSfBp')
		json_buck = conn.get_bucket('rep-seq-jsons')
		json_names = [str(key.name) for key in json_buck.list()]

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
					#self.SHM_DF = [key, group['name'], sample, 'error']
					#self.SHM_DF['Data'].loc[sample] = SHM_data

		self.SHM_DF = pd.DataFrame(ls, columns = SHM_columns)

		#Creates a Panda DataFrame for the clone information. The columns are [Group  Top 1  Top 10  Top 100]
		Clone_columns = ['Group', 'Top 1', 'Top 10', 'Top 100']

		self.clone_DF = pd.DataFrame(columns=Clone_columns)

		for igroup, group in enumerate(self.groups):
		    for sample in group['samples']:
		        top = self.groups[igroup]['sample data'][sample]['top_clone_fraction']
		        top_10 = self.groups[igroup]['sample data'][sample]['top_10_fraction']
		        top_100 = self.groups[igroup]['sample data'][sample]['top_100_fraction']
		        
		        self.clone_DF.loc[sample] = [group['name'], top, top_10, top_100]
		        

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

		self.clone_DF.plot(kind='bar', stacked='True');

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
