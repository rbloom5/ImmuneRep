import Experiment_class
reload(Experiment_class)
from Experiment_class import Experiment
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sets
from matplotlib_venn import venn2
import numpy as np
import math

	

def christianDF():
	datain = np.loadtxt('patient_index.txt', dtype='string', delimiter='/t')

	columns = ['Patient', 'Diagnosis', 'Sex', 'Age', 'EDSS', 'Relapse', 'Source', 'File']

	ls = []
	for i,row in enumerate(datain):
	    patient = math.ceil(i/2)+1
	    filenumber =  row.split('\t')[0]
	    source = row.split('\t')[1].split('-')[1]
	    diagnosis = row.split('\t')[1].split('-')[0]
	    sex = row.split('\t')[6].split(',')[0]
	    age = row.split('\t')[6].split(',')[1].split(' ')[2]
	    EDSS =  row.split('\t')[6].split(',')[2].split('=')[1]
	    relapse = row.split('\t')[6].split(',')[3]
	    
	    ls.append([patient, diagnosis, sex, age, EDSS, relapse, source, filenumber])
	    
	return pd.DataFrame(ls, columns=columns)

def create_linked_list():

	DF = christianDF()

	ls = []
	for index, row in sDF.iterrows():
	    ls.append([str(row[0])+" "+row[6], row[7]])
	print ls

	data = [ls[0]+ls[1],ls[2]+ls[3],ls[4]+ls[5],ls[7]+ls[6],ls[9]+ls[8],ls[11]+ls[10],ls[13]+ls[12],ls[15]+ls[14]]

	linked = {}
	for patient in data:
	    group_list = []
	    
	    dict1 = {}
	    dict1["name"] = patient[0]
	    dict1["samples"] = [patient[1]]
	    
	    dict2 = {}
	    dict2["name"] = patient[2]
	    dict2["samples"] = [patient[3]]
	    
	    group_list.append(dict1)
	    group_list.append(dict2)
	    
	    linked[patient[0].split(' ')[0].split('.')[0]] = Experiment(group_list)

	return linked