#!/usr/bin/python

import vdj_fasta
from vdj_fasta import multiprocess_vdj

import rep_seq
from rep_seq import Rep_sequence_analysis

import subprocess
import time
import pickle
import os


##################################################
## Main Pipeline for processing immune rep data ##
# must be run on amazon ubuntu machine with all the dependencies built in (see wiki)

## 1.)  input names of fasta files from patient_repertoire _data
## folder in s3 (do not put .fasta extension) as list

## 2.)  Set options

## the script will do the rest of the work


def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output




#############################
### SET OPTIONS AND INPUTS ##
#############################

fileIDs = [	#'SRR1383453',\
			#'SRR1383472',\
			'SRR1383450',\
			'SRR1383455',\
			'SRR1383466',\
			'SRR1383470',\
			'SRR1383476',\
			] #must be a list of file id's (no .fasta extension)


num_sequences = 100000 # max number of sequences to analyze. 50 - 100k is a good depth
ABtype = None # set if you only want to analyze a certain class of antibody (prob not)
plots = False #if you are running on a local machine and want to see plots - set to True




###########################
####### MAIN SCRIPT #######
###########################


# Run VDJ-fasta - store output in s3 folder: ''
#run_vdjfasta(fileIDs, num_sequences)

# load in vdj-fasta output and pull out features using rep-seq
ext = '.VDJ.H3.L3.CH1.fa'
s3_dir = 's3://vdjfasta-output/'
counter = 0

for f in fileIDs:
	start = time.time()
	print "started processing", f
	f_string = s3_dir + f + ext

	#copy from s3 to local
	bash('aws s3 cp %s /home/ubuntu/tempvdj/%s'%(f_string,f+ext))

	# use rep seq to get all properties
	Rep = Rep_sequence_analysis.Rep_seq(['/home/ubuntu/tempvdj/'+f+ext],num=num_sequences)
	Rep.find_clones(parallel=True)
	print "calculating statistics..."
	Rep.get_stats()
	Rep.tree()
	if plots:
		Rep.plots()

	#output features and Rep-seq object to files
	print "saving features to json..."
	Rep.output_features('/home/ubuntu/tempvdj/'+f+'.json')
	with open('/home/ubuntu/tempvdj/'+f+'.pkl', 'wb') as output:
		pickle.dump(Rep, output, pickle.HIGHEST_PROTOCOL)

	#Copy all files to the appropriate folder in s3
	print "copying files to s3"
	bash('aws s3 cp /home/ubuntu/tempvdj/'+f+'.json s3://rep-seq-jsons/')
	bash('aws s3 cp /home/ubuntu/tempvdj/'+f+'.pkl s3://rep-seq-objects/')

	
	node_files = [i for i in os.listdir('/home/ubuntu/tree_output') if i.endswith('.fasta.node.txt')]
	for nf in node_files:
		aws_loc = 's3://tree-output/'+f+'/' + nf
		bash('aws s3 cp /home/ubuntu/tree_output/'+nf + ' '+aws_loc)

	#remove all temp files that have been transferred
	os.system('sudo rm -f -r /home/ubuntu/tempvdj/*')
	os.system('sudo rm -f -r /home/ubuntu/tree_output/*')
	os.system('sudo rm -f -r /home/ubuntu/parsed_fasta/*')

	counter+=1
	print f, "complete!  %s of %s files finished\n"%(counter,len(fileIDs))
	print "processing time:", time.time()-start










