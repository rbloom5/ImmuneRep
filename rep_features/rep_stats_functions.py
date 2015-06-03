#!/usr/bin/python
import numpy as np
from ete2 import Tree
import pandas as pd



def Make_VJ_Matrix():
    #format for calling location matrix_name['J_name'].loc['V_name']

<<<<<<< Updated upstream
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



    J_names = ['IGHJ1', 'IGHJ1P', 'IGHJ2', 'IGHJ2P', 'IGHJ3', 'IGHJ3', 'IGHJ3P', 'IGHJ3P', 'IGHJ4', 'IGHJ4', 'IGHJ4', 'IGHJ5', 'IGHJ5', 'IGHJ6',\
               'IGHJ6', 'IGHJ6', 'IGHJ6', 'mIGHJ1', 'mIGHJ1', 'mIGHJ1', 'mIGHJ2', 'mIGHJ2', 'mIGHJ3', 'mIGHJ3', 'mIGHJ4', 'mIGKJ1', 'mIGKJ1', 'mIGKJ2', \
               'mIGKJ2', 'mIGKJ2', 'mIGKJ3', 'mIGKJ3', 'mIGKJ4', 'mIGKJ4', 'mIGKJ5', 'mIGLJ1', 'mIGLJ2', 'mIGLJ3', 'mIGLJ3P', 'mIGLJ4', 'IGLJ1', 'IGLJ2/3', \
               'IGLJ2/3', 'IGLJ2/3', 'IGLJ4', 'IGLJ5', 'IGLJ5', 'IGLJ6', 'IGLJ7', 'IGLJ7', 'IGKJ1', 'IGKJ2', 'IGKJ2', 'IGKJ2', 'IGKJ2', 'IGKJ3', 'IGKJ4', \
               'IGKJ4', 'IGKJ5', 'TRAJ1', 'TRAJ10', 'TRAJ11', 'TRAJ12', 'TRAJ13', 'TRAJ13', 'TRAJ14', 'TRAJ15', 'TRAJ15', 'TRAJ16', 'TRAJ17', 'TRAJ18', 'TRAJ19', \
               'TRAJ2', 'TRAJ20', 'TRAJ21', 'TRAJ22', 'TRAJ23', 'TRAJ23', 'TRAJ24', 'TRAJ24', 'TRAJ25', 'TRAJ26', 'TRAJ27', 'TRAJ28', 'TRAJ29', 'TRAJ3', 'TRAJ30', \
               'TRAJ31', 'TRAJ32', 'TRAJ32', 'TRAJ33', 'TRAJ34', 'TRAJ35', 'TRAJ36', 'TRAJ37', 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRAJ4', 'TRAJ40', 'TRAJ41', 'TRAJ42', \
               'TRAJ43', 'TRAJ44', 'TRAJ45', 'TRAJ46', 'TRAJ47', 'TRAJ47', 'TRAJ48', 'TRAJ49', 'TRAJ5', 'TRAJ50', 'TRAJ51', 'TRAJ52', 'TRAJ53', 'TRAJ54', 'TRAJ55', \
               'TRAJ56', 'TRAJ57', 'TRAJ58', 'TRAJ59', 'TRAJ6', 'TRAJ60', 'TRAJ61', 'TRAJ7', 'TRAJ8', 'TRAJ9', 'TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', \
               'TRBJ1-6', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-2P', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7', 'TRBJ2-7', 'TRGJ1', 'TRGJ1', 'TRGJ2', 'TRGJP', \
               'TRGJP1', 'TRGJP2', 'TRDJ1', 'TRDJ2', 'TRDJ3', 'TRDJ4']
=======
    V_names = ['IGHV1-12', 'IGHV4-4', 'IGHV1-17', 'IGHV1-14', 'IGHV3-15', 'IGHV1-18', 'IGHV3-11', 'IGHV3-13', 'IGHV(III)-16-1', \
                'IGHV3-33-2', 'IGHV4-30-4%31', 'IGHV3-19', 'IGHV(II)-33-1', 'IGHV3-7', 'IGHV3-6', 'IGHV3-30-2', 'IGHV(III)-5-2', \
                'IGHV(II)-30-1', 'IGHV4-59%61', 'IGHV(III)-25-1', 'IGHV(IV)-44-1', 'IGHV(II)-65-1', 'IGHV1-24', 'IGHV3-16', 'IGHV(III)-2-1', \
                'IGHV1-NL1', 'IGHV2-70', 'IGHV4-30-2', 'IGHV3-63', 'IGHV3-64', 'IGHV3-65', 'IGHV4-80', 'IGHV7-4-1', 'IGHV(III)-47-1', \
                'IGHV5-78', 'IGHV3-d', 'IGHV(II)-74-1', 'IGHV3-h', 'IGHV(III)-11-1', 'IGHV3-79', 'IGHV3-73', 'IGHV3-72', 'IGHV3-71', \
                'IGHV3-76', 'IGHV3-75', 'IGHV3-74', 'IGHV(III)-82', 'IGHV1-8', 'IGHV7-34-1%40%NL1.rn', 'IGHV1-3', 'IGHV1-2', 'IGHV4-39', \
                'IGHV4-34', 'IGHV7-81', 'IGHV(II)-62-1', 'IGHV(II)-51-2', 'IGHV3-48', 'IGHV3-49', 'IGHV3-47', 'IGHV3-42', 'IGHV(II)-22-1', \
                'IGHV3-41', 'IGHV3-60%62', 'IGHV7-27', 'IGHV3-9', 'IGHV(II)-44-2', 'IGHV1-f', 'IGHV(III)-38-1', 'IGHV1-c', 'IGHV6-1', \
                'IGHV3-52', 'IGHV3-54', 'IGHV3-57', 'IGHV(III)-67-4', 'IGHV(III)-67-3', 'IGHV(II)-31-1', 'IGHV2-5', 'IGHV4-55', 'IGHV3-43', \
                'IGHV(II)-53-1', 'IGHV1-45', 'IGHV1-46', 'IGHV3-25', 'IGHV3-20', 'IGHV3-21', 'IGHV3-22', 'IGHV3-23', 'IGHV3-29', 'IGHV(III)-51-1', \
                'IGHV(II)-78-1', 'IGHV1-58', 'IGHV5-a', 'IGHV5-51', 'IGHV(II)-15-1', 'IGHV3-37', 'IGHV3-36', 'IGHV3-35', 'IGHV3-32', 'IGHV2-26', \
                'IGHV3-38', 'IGHV(III)-13-1', 'IGHV1-67', 'IGHV(III)-44', 'IGHV1-68', 'IGHV1-69', 'IGHV(II)-28-1', 'IGHV(II)-60-1', 'IGHV(II)-49-1', \
                'IGHV(III)-26-1', 'IGHV7-56', 'IGHV4-28', 'IGHV3-53%66', 'IGHV(II)-20-1', 'IGHV2-10', 'IGHV3-30%33rn', 'IGHV4-b', 'IGHV(III)-76-1', \
                'mIGHV6-3%6', 'IGHV3-NL1', 'mIGHV7-3']


    J_names = ['TRAJ57', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7', 'TRBJ2-1', 'TRBJ2-2', 'IGHJ1P', 'TRAJ30', 'TRAJ31', 'TRAJ32', 'TRAJ33', \
               'TRAJ34', 'TRAJ35', 'TRAJ36', 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRBJ2-2P', 'IGKJ1', 'IGKJ2', 'IGKJ3', 'IGKJ4', 'IGKJ5', 'TRBJ1-6', \
               'TRBJ1-5', 'TRBJ1-4', 'TRBJ1-3', 'TRBJ1-2', 'TRBJ1-1', 'TRAJ29', 'TRAJ28', 'TRAJ27', 'TRAJ26', 'TRAJ25', 'TRAJ24', 'TRAJ23', 'TRAJ22', \
               'TRAJ21', 'TRAJ20', 'mIGKJ4', 'mIGKJ5', 'mIGKJ1', 'mIGKJ2', 'mIGKJ3', 'mIGLJ3P', 'TRBJ2-3', 'TRAJ58', 'TRAJ59', 'TRAJ52', 'TRAJ53', \
               'TRAJ50', 'TRAJ51', 'TRAJ56', 'IGLJ2/3', 'TRAJ54', 'TRAJ55', 'TRGJP', 'TRDJ4', 'TRDJ1', 'TRDJ3', 'TRDJ2', 'TRAJ49', 'TRAJ48', 'TRAJ41', \
               'TRAJ40', 'TRAJ43', 'TRAJ42', 'TRAJ45', 'TRAJ44', 'TRAJ47', 'TRAJ46', 'IGLJ1', 'IGLJ7', 'IGLJ6', 'IGLJ5', 'IGLJ4', 'mIGHJ3', 'mIGHJ2',\
               'mIGHJ1', 'mIGHJ4', 'IGHJ2P', 'TRGJ1', 'TRGJ2', 'IGHJ6', 'IGHJ5', 'IGHJ4', 'IGHJ3', 'IGHJ2', 'IGHJ1', 'TRAJ4', 'TRAJ5', 'TRAJ6', 'TRAJ7', \
               'TRAJ1', 'TRAJ2', 'TRAJ3', 'TRAJ8', 'TRAJ9', 'TRAJ61', 'TRAJ60', 'mIGLJ4', 'mIGLJ3', 'mIGLJ2', 'mIGLJ1', 'TRAJ16', 'TRAJ17', 'TRAJ14',\
               'TRAJ15', 'TRAJ12', 'TRAJ13', 'TRAJ10', 'TRAJ11', 'TRAJ18', 'TRAJ19', 'IGHJ3P', 'TRGJP2', 'TRGJP1']
>>>>>>> Stashed changes


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

def avg_node_mutation_length(tree_dict):
    
    tree_mutation_dict = {}
    
    for key in tree_dict:
        ls = []
        for node in tree_dict[key].iter_descendants():
            ls.append(len(node.mutations))
        tree_mutation_dict[key] = [np.mean(ls), np.std(ls)]
        
    return tree_mutation_dict

def generations(tree_dict):
    node_generations = {}
    for key in tree_dict:
        node_generations[key] = tree_dict[key].get_farthest_node()[1]
    return node_generations


def calculate_vj_shm(vj_clones):
    sh_means = {}
    sh_stdevs = {}
    for vj in vj_clones:
        # print vj_clones[vj]
        all_sh = []
        for clone in vj_clones[vj]:
            all_sh.append(vj_clones[vj][clone].sh)
        sh_means[vj] = np.array(all_sh).mean() 
        sh_stdevs[vj] = np.array(all_sh).std() 

    return sh_means, sh_stdevs










