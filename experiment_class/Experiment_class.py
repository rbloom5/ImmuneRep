import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 

class Experiment:
	def __init__(self,groups,datadir='',defaultgroupspath='defaultgroups.txt'):
		defaultgroups=json.load(open(defaultgroupspath)) 

		#Checks for groups with no input samples and looks for a default
		for igroup, group in enumerate(groups):
			if 'samples' not in group:
				groups[igroup]['samples']=defaultgroups[group['name']]

		#Pulls in the sample data for each sample
		for igroup, group in enumerate(groups):
			groups[igroup]['sample data'] = {}
			for sample in group['samples']:
				groups[igroup]['sample data'][sample]=json.load(open(datadir+'/'+sample+'.txt'))

		self.groups = groups

		#Creates a color dictionary that fades per sample and is based on the input color
		self.color_dict = {}
		for group in self.groups:
		    color = sns.color_palette(group['color'], len(group['samples']))
		    for i, sample in enumerate(group['samples']):
		        self.color_dict[sample] = color[i]

		#Creates two Pandas DataFrames with the columsn [Antibody  Group  Sample  Data]. The first repersents the data column as a list of
		#somatic hypermutations. The second creates a new row for every SHM.
		SHM_columns = ['Antibody', 'Group', 'Sample', 'Data']
		antibody_class_list = ['all classes', 'IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD']

		ls_split = []
		ls = []
		for key in antibody_class_list:
			for igroup,group in enumerate(self.groups):
				for sample in group['samples']:
					ls.append([key, group['name'], sample, self.groups[igroup]['sample data'][sample]['sh_dict'][key]])
					for i in self.groups[igroup]['sample data'][sample]['sh_dict'][key]:
						ls_split.append([key, group['name'], sample, i])
		DF_split = pd.DataFrame(ls_split, columns = SHM_columns)
		DF = pd.DataFrame(ls, columns = SHM_columns)

		self.SHM_DF = DF
		self.SHM_DF_split = DF_split

		#Creates a Panda DataFrame for the clone information. The columns are [Group  Top 1  Top 10  Top 100]
		Clone_columns = ['Group', 'Top 1', 'Top 10', 'Top 100']

		self.Clone_DF = pd.DataFrame(columns=Clone_columns)

		for igroup, group in enumerate(self.groups):
		    for sample in group['samples']:
		        top = self.groups[igroup]['sample data'][sample]['top_clone_fraction']
		        top_10 = self.groups[igroup]['sample data'][sample]['top_10_fraction']
		        top_100 = self.groups[igroup]['sample data'][sample]['top_100_fraction']
		        
		        self.Clone_DF.loc[sample] = [group['name'], top, top_10, top_100]
		        

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

		self.Clone_DF.plot(kind='bar', stacked='True');
