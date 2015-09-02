#!/usr/bin/python

import vdj_fasta
from vdj_fasta import multiprocess_vdj_local

import rep_seq
from rep_seq import Rep_sequence_analysis

import subprocess
import time
import pickle
import os
from Bio import SeqIO


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


def make_dir(dirname):
	if not os.path.exists('%s/'%dirname):
		os.makedirs('%s/'%dirname)
	else:
		os.system('rm %s/*'%dirname)

def copy_to_s3(files, local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://%s/'%s3_dir + f
		bash('aws s3 cp %s/'%local_dir+f+' '+aws_loc)

def copy_from_s3(files,local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://%s/'%s3_dir + f
		bash('aws s3 cp ' + aws_loc + ' %s/'%local_dir + f)

def fastq_dump(f, clean=True):
	os.system('./packages/sra_toolkit/bin/fastq-dump %s'%f)
	copy_to_s3(f+'.fastq', '.', 'patient_repertoire_data')
	if clean:
		SeqIO.convert(f+'.fastq', 'fastq', f+'.fasta', 'fasta')
		copy_to_s3(f+'.fasta', '.', 'clean-repertoire-data')



#############################
### SET OPTIONS AND INPUTS ##
#############################

fileIDs = [	#'SRR1383453',\
			

			# # 'newQuake_MID0-len-q-t',\
			# # 'newQuake_MID1-len-q-t',\
			# # 'newQuake_MID2-len-q-t',\
			# # 'newQuake_MID6-len-q-t',\
			# # 'newQuake_MID7-len-q-t',\
			# # 'newQuake_MID8-len-q-t',\
			# # 'newQuake_MID9-len-q-t',\
			# # 'newQuake_post_vax_MID0-len-q-t',\
			# # 'newQuake_post_vax_MID1-len-q-t',\
			# 'newQuake_post_vax_MID2-len-q-t',\
			# 'newQuake_post_vax_MID6-len-q-t',\
			# 'newQuake_post_vax_MID7-len-q-t',\
			# 'newQuake_post_vax_MID8-len-q-t',\
			# 'newQuake_post_vax_MID9-len-q-t',\
			# 'SRR1298742',\
			# 'SRR1383448',\
			# 'SRR1383453',\
			# 'SRR1383463',\
			# 'SRR1383472',\
			# 'SRR1383450',\
			# 'SRR1383455',\
			# 'SRR1383466',\
			# 'SRR1383476',\

			# #christian
			# 'SRR1298740',\
			# 'SRR1297001',\

			# #allergies
			# 'SRR1171336',\
			# 'SRR1171337',\
			# 'SRR1171338',\
			# 'SRR1171339',\
			# 'SRR1171340',\
			# 'SRR1171341',\
			# 'SRR1171342',\
			# 'SRR1171343',\
			# # 'SRR1171344',\
			# 'SRR1171345',\


			# #b_cell subsets
			# 'SRR1168779',\
			# 'SRR1168788',\
			# 'SRR1168789',\
			# 'SRR1168790',\
			# 'SRR1168792',\
			# 'SRR1168794',\

			#more christian
			# 'SRR1298383',\
			# 'SRR1298730',\
			# 'SRR1298731',\
			# 'SRR1298732',\
			# 'SRR1298733',\
			# 'SRR1298734',\
			# 'SRR1298735',\
			# 'SRR1298736',\
			# 'SRR1298737',\
			# 'SRR1298738',\
			# 'SRR1298739',\
			# 'SRR1298741',\
			# 'SRR1298743',\

			#SLE
			# 'SRR1959703',\
			# 'SRR1960371',\
			# 'SRR1961400',\
			# 'SRR1964710',\
			# 'SRR1964711',\
			# 'SRR1964712',\
			'SRR1964713',\
			'SRR1964786',\
			'SRR1964787',\
			'SRR1964788',\
			'SRR1964792',\
			'SRR1964793',\
			'SRR1964794',\
			'SRR1964795',\
			'SRR1964796',\
			'SRR1964797',\
			'SRR1964798',\
			'SRR1964799',\
			'SRR1964800',\
			'SRR1964801',\


			] #must be a list of file id's (no .fasta extension)


num_sequences = 100000 # max number of sequences to analyze. 50 - 100k is a good depth
ABtype = None # set if you only want to analyze a certain class of antibody (prob not)
plots = False #if you are running on a local machine and want to see plots - set to True




###########################
####### MAIN SCRIPT #######
###########################

# #get files, if necessary
# for ids in fileIDs:
# 	print "downloading %s from SRA"%ids
# 	counter=1
# 	# os.system('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra'%(ids[:6], ids, ids))
# 	fastq_dump(ids)
# 	os.system('rm -r ncbi')
# 	os.system('rm %s.fastq'%ids)
# 	os.system('rm %s.fasta'%ids)


# # Run VDJ-fasta - store output in s3 folder: ''
# multiprocess_vdj_local.run_vdjfasta(fileIDs, num_sequences)

# load in vdj-fasta output and pull out features using rep-seq
ext = '.VDJ.H3.L3.CH1.fa'
s3_dir = 's3://vdjfasta-output/'
counter = 0

for f in fileIDs:
	start = time.time()
	print "started processing", f
	f_string = s3_dir + f + ext

	make_dir('/home/ubuntu/tempvdj')


	#copy from s3 to local
	bash('aws s3 cp %s /home/ubuntu/tempvdj/%s'%(f_string,f+ext))
	# bash('aws s3 cp s3://tree-output/'+f+'/' ' /home/ubuntu/tree_output/ --recursive')


	# # use rep seq to get all properties
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
	os.system('sudo rm -f -r /home/ubuntu/tree_output')
	os.system('sudo rm -f -r /home/ubuntu/parsed_fasta')



	#remove all temp files that have been transferred
	os.system('sudo rm -f -r /home/ubuntu/tempvdj/*')


	counter+=1
	print f, "complete!  %s of %s files finished\n"%(counter,len(fileIDs))
	print "processing time:", time.time()-start

os.system('sudo rm -f -r /home/ubuntu/tempvdj')











