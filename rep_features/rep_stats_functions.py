#!/usr/bin/python
import numpy
from ete2 import Tree
from pandas import *



def Make_VJ_Matrix():
    #format for calling location matrix_name['J_name'].loc['V_name']

    V_names = ['IGHV1-2', 'IGHV1-3', 'IGHV1-8', 'IGHV1-12', 'IGHV1-14', 'IGHV1-17', 'IGHV1-18', 'IGHV1-24', 'IGHV1-38', 'IGHV1-45', 'IGHV1-46',\
               'IGHV1-58', 'IGHV1-67', 'IGHV1-68', 'IGHV1-69', 'IGHV2-5', 'IGHV2-10', 'IGHV2-26', 'IGHV2-70', 'IGHV3-6', 'IGHV3-7', 'IGHV3-9', \
               'IGHV3-11', 'IGHV3-13', 'IGHV3-15', 'IGHV3-16', 'IGHV3-19', 'IGHV3-20', 'IGHV3-21', 'IGHV3-22', 'IGHV3-23', 'IGHV3-25', 'IGHV3-29', \
               'IGHV3-30', 'IGHV3-32', 'IGHV3-33', 'IGHV3-35', 'IGHV3-36', 'IGHV3-37', 'IGHV3-38', 'IGHV3-41', 'IGHV3-42', 'IGHV3-43', 'IGHV3-47', \
               'IGHV3-48', 'IGHV3-49', 'IGHV3-50', 'IGHV3-52', 'IGHV3-53', 'IGHV3-54', 'IGHV3-57', 'IGHV3-60', 'IGHV3-62', 'IGHV3-63', 'IGHV3-64', \
               'IGHV3-65', 'IGHV3-66', 'IGHV3-69', 'IGHV3-71', 'IGHV3-72', 'IGHV3-73', 'IGHV3-74', 'IGHV3-75', 'IGHV3-76', 'IGHV3-79', 'IGHV4-4', \
               'IGHV4-28', 'IGHV4-30', 'IGHV4-31', 'IGHV4-34', 'IGHV4-38', 'IGHV4-39', 'IGHV4-55', 'IGHV4-59', 'IGHV4-61', 'IGHV4-80', 'IGHV5-10', \
               'IGHV5-51', 'IGHV5-78', 'IGHV6-1', 'IGHV7-4', 'IGHV7-27', 'IGHV7-34', 'IGHV7-40', 'IGHV7-56', 'IGHV7-77', 'IGHV7-81']

    J_names = ['IGHJ1','IGHJ2','IGHJ3','IGHJ4','IGHJ5','IGHJ6']

    return DataFrame(np.zeros((len(V_names),len(J_names))),index=V_names,columns=J_names)

def calculate_tree_size(rep_obj):

	tree_size_DF = Make_VJ_Matrix()

	for vj_pair in rep_obj.tree_dict:
		V = vj_pair.split('_').[0]
		J = vj_pair.split('_').[1]

		tree_size_DF[J].loc[V] = len(rep_obj.tree_dict[vj_pair].get_descendants())

	return tree_size_DF


def calculate_something(Rep):

	#calculate something useful on the rep

	return 'the useful thing'


def calculate_something_else(Rep):

	# calculate something else useful

	return 'the other useful thing'