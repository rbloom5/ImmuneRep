#!/usr/bin/python

import vdj_fasta
import rep_seq
import subprocess
import time
import pickle


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

fileIDs = [] #must be a list of file id's (no .fasta extension)
num_sequences = 100000 # max number of sequences to analyze. 50 - 100k is a good depth
ABtype = None # set if you only want to analyze a certain class of antibody (prob not)
plots = False #if you are running on a local machine and want to see plots - set to True




###########################
####### MAIN SCRIPT #######
###########################


# Run VDJ-fasta - store output in s3 folder: ''
run_vdjfasta(fileIDs)

# load in vdj-fasta output and pull out features using rep-seq
ext = '.VDJ.H3.L3.CH1.fa'
s3_dir = 's3://vdjfasta-output/'
counter = 0
for f in fileIDs:
	print "started processing", f
	f_string = s3_dir + f + ext

	#copy from s3 to local
	bash('aws s3 cp %s /home/ubuntu/tempvdj'%f_string)

	# use rep seq to get all properties
	Rep = Rep_seq(['/home/ubuntu/tempvdj/'+f+ext],num=50000)
	Rep.find_clones(parallel=True)
	Rep.get_stats()
	Rep.tree()
	if plots:
		Rep.plots()

	#output features and Rep-seq object to files
	Rep.output_features('/home/ubuntu/tempvdj/'+f+'.json')
	with open('home/ubuntu/tempvdj/'+f+'.pkl', 'wb') as output:
		pickle.dump(Rep, output, pickle.HIGHEST_PROTOCOL)

	#Copy all files to the appropriate folder in s3
	bash('aws s3 cp /home/ubuntu/tempvdj/'+f+'.json s3://rep-seq-jsons/')
	bash('aws s3 cp /home/ubuntu/tempvdj/'+f+'.pkl s3://rep-seq-objects/')

	bash('aws s3 mkdir s3://tree-output/'+f+'/')
	bash('aws s3 cp /home/ubuntu/tree_output/*.fasta.node.txt'+f+'.pkl s3://tree-output/'+f+'/')

	#remove all temp files that have been transferred
	bash('rm -f /home/ubuntu/tempvdj/*')
	bash('rm -f /home/ubuntu/tree_output/*')
	bash('rm -f /home/ubuntu/parsed_fasta/*')

	counter+=1
	print f, "complete!  %s of %s files finished\n"%(counter,len(fileIDs))










