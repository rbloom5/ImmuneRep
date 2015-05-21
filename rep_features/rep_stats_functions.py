#!/usr/bin/python
import numpy as np
from ete2 import Tree
import pandas as pd



def Make_VJ_Matrix():
    #format for calling location matrix_name['J_name'].loc['V_name']

    V_names = ['IGHV(II)-15-1', 'IGHV(II)-20-1', 'IGHV(II)-22-1', 'IGHV(II)-28-1', 'IGHV(II)-30-1', 'IGHV(II)-31-1', 'IGHV(II)-33-1', \
    			'IGHV(II)-44-2', 'IGHV(II)-49-1', 'IGHV(II)-51-2', 'IGHV(II)-53-1', 'IGHV(II)-60-1', 'IGHV(II)-62-1', 'IGHV(II)-65-1', \
    			'IGHV(II)-74-1', 'IGHV(II)-78-1', 'IGHV(III)-11-1', 'IGHV(III)-13-1', 'IGHV(III)-16-1', 'IGHV(III)-2-1', 'IGHV(III)-25-1', \
    			'IGHV(III)-26-1', 'IGHV(III)-38-1', 'IGHV(III)-44', 'IGHV(III)-47-1', 'IGHV(III)-5-2', 'IGHV(III)-51-1', 'IGHV(III)-67-3', \
    			'IGHV(III)-67-4', 'IGHV(III)-76-1', 'IGHV(III)-82', 'IGHV(IV)-44-1', 'IGHV1-12', 'IGHV1-14', 'IGHV1-17', 'IGHV1-17', \
    			'IGHV1-18', 'IGHV1-18', 'IGHV1-2', 'IGHV1-2', 'IGHV1-2', 'IGHV1-2', 'IGHV1-24', 'IGHV1-3', 'IGHV1-3', 'IGHV1-45', 'IGHV1-45',\
    			'IGHV1-45', 'IGHV1-46', 'IGHV1-46', 'IGHV1-46', 'IGHV1-58', 'IGHV1-58', 'IGHV1-67', 'IGHV1-67', 'IGHV1-68', 'IGHV1-69', \
    			'IGHV1-69', 'IGHV1-69', 'IGHV1-69', 'IGHV1-69', 'IGHV1-69', 'IGHV1-69', 'IGHV1-69', 'IGHV1-69', 'IGHV1-69', 'IGHV1-69', \
    			'IGHV1-69', 'IGHV1-69', 'IGHV1-8', 'IGHV1-NL1', 'IGHV1-c', 'IGHV1-f', 'IGHV1-f', 'IGHV2-10', 'IGHV2-26', 'IGHV2-5', 'IGHV2-5',\
    			'IGHV2-5', 'IGHV2-5', 'IGHV2-5', 'IGHV2-5', 'IGHV2-5', 'IGHV2-5', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70',\
    			'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV2-70', 'IGHV3-11', 'IGHV3-11', \
    			'IGHV3-11', 'IGHV3-13', 'IGHV3-13', 'IGHV3-13', 'IGHV3-15', 'IGHV3-15', 'IGHV3-15', 'IGHV3-15', 'IGHV3-15', 'IGHV3-15',\
    			'IGHV3-15', 'IGHV3-15', 'IGHV3-16', 'IGHV3-16', 'IGHV3-19', 'IGHV3-20', 'IGHV3-21', 'IGHV3-21', 'IGHV3-22', 'IGHV3-22', \
    			'IGHV3-23', 'IGHV3-23', 'IGHV3-23', 'IGHV3-23', 'IGHV3-23', 'IGHV3-25', 'IGHV3-25', 'IGHV3-25', 'IGHV3-29', 'IGHV3-30$33rn',\
    			'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', \
    			'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', \
    			'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30-2', 'IGHV3-30$33rn', 'IGHV3-30$33rn', \
    			'IGHV3-32', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-30$33rn', 'IGHV3-33-2', 'IGHV3-35', \
    			'IGHV3-36', 'IGHV3-36', 'IGHV3-37', 'IGHV3-37', 'IGHV3-38', 'IGHV3-38', 'IGHV3-41', 'IGHV3-42', 'IGHV3-42', 'IGHV3-42', \
    			'IGHV3-43', 'IGHV3-43', 'IGHV3-47', 'IGHV3-47', 'IGHV3-47', 'IGHV3-48', 'IGHV3-48', 'IGHV3-48', 'IGHV3-49', 'IGHV3-49',\
    			'IGHV3-49', 'IGHV3-49', 'IGHV3-49', 'IGHV3-52', 'IGHV3-52', 'IGHV3-52', 'IGHV3-53$66', 'IGHV3-53$66', 'IGHV3-53$66', \
    			'IGHV3-54', 'IGHV3-54', 'IGHV3-54', 'IGHV3-54', 'IGHV3-57', 'IGHV3-57', 'IGHV3-6', 'IGHV3-60$62', 'IGHV3-60$62', \
    			'IGHV3-63', 'IGHV3-63', 'IGHV3-64', 'IGHV3-64', 'IGHV3-64', 'IGHV3-64', 'IGHV3-64', 'IGHV3-65', 'IGHV3-65', 'IGHV3-65', \
    			'IGHV3-53$66', 'IGHV3-53$66', 'IGHV3-53$66', 'IGHV3-53$66', 'IGHV3-7', 'IGHV3-7', 'IGHV3-71', 'IGHV3-72', 'IGHV3-73', \
    			'IGHV3-73', 'IGHV3-74', 'IGHV3-74', 'IGHV3-74', 'IGHV3-75', 'IGHV3-76', 'IGHV3-76', 'IGHV3-79', 'IGHV3-9', 'IGHV3-d', \
    			'IGHV3-h', 'IGHV3-h', 'IGHV4-28', 'IGHV4-28', 'IGHV4-28', 'IGHV4-28', 'IGHV4-28', 'IGHV4-30-2', 'IGHV4-30-2', 'IGHV4-30-2', \
    			'IGHV4-30-2', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', \
    			'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', \
    			'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-30-4$31', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', \
    			'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-34', 'IGHV4-39', 'IGHV4-39', 'IGHV4-39', \
    			'IGHV4-39', 'IGHV4-39', 'IGHV4-39', 'IGHV4-39', 'IGHV4-4', 'IGHV4-4', 'IGHV4-4', 'IGHV4-4', 'IGHV4-4', 'IGHV4-4', 'IGHV4-4', \
    			'IGHV4-55', 'IGHV4-55', 'IGHV4-55', 'IGHV4-55', 'IGHV4-55', 'IGHV4-55', 'IGHV4-55', 'IGHV4-55', 'IGHV4-55', 'IGHV4-59$61', \
    			'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', \
    			'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', 'IGHV4-59$61', \
    			'IGHV4-59$61', 'IGHV4-80', 'IGHV4-b', 'IGHV4-b', 'IGHV5-51', 'IGHV5-51', 'IGHV5-51', 'IGHV5-51', 'IGHV5-51', 'IGHV5-78', \
    			'IGHV5-78', 'IGHV5-a', 'IGHV5-a', 'IGHV5-a', 'IGHV5-a', 'IGHV6-1', 'IGHV6-1', 'IGHV7-27', 'IGHV7-34-1$40$NL1.rn', 'IGHV7-4-1', \
    			'IGHV7-4-1', 'IGHV7-4-1', 'IGHV7-34-1$40$NL1.rn', 'IGHV7-34-1$40$NL1.rn', 'IGHV7-56', 'IGHV7-81', 'IGHV7-34-1$40$NL1.rn', \
    			'IGHV7-34-1$40$NL1.rn', 'IGHV7-34-1$40$NL1.rn', 'IGHV7-34-1$40$NL1.rn', 'IGHV7-34-1$40$NL1.rn']


    J_names = ['IGHJ1','IGHJ2','IGHJ3','IGHJ4','IGHJ5','IGHJ6']

    return pd.DataFrame(np.zeros((len(V_names),len(J_names))),index=V_names,columns=J_names)

def calculate_tree_size(rep_obj, pruned=False):

	if pruned==False: dictionary = rep_obj.tree_dict
	if pruned==True: dictionary = rep_obj.pruned_tree_dict

	tree_size_DF = Make_VJ_Matrix()

	for vj_pair in rep_obj.tree_dict:
		V = vj_pair.split('_')[0]
		J = vj_pair.split('_')[1]


		tree_size_DF[J].loc[V] = len(dictionary[vj_pair].get_descendants())

	tree_size_dict = tree_size_DF.to_dict()

	return tree_size_dict


def calculate_cdr_lengths(obj_dict, cdr='cdr3'):
	len_dict = {}
	for obj in obj_dict:
		cdr_len = len(str(obj_dict[obj].cdr3))
		if cdr_len not in len_dict:
			len_dict[cdr_len] = 1
		else:
			len_dict[cdr_len]+=1
					
	return len_dict







