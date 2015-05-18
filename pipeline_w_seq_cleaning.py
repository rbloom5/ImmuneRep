#!/usr/bin/python

import vdj_fasta
from vdj_fasta import multiprocess_vdj_local

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

#############################
### SET OPTIONS AND INPUTS ##
#############################

fileIDs = [	\
			'SRR747232-58-60',\
			'SRR747766',\
			'SRR747785',\

			] #must be a list of file id's (no .fasta extension)


num_sequences = 100000 # max number of sequences to analyze. 50 - 100k is a good depth
ABtype = None # set if you only want to analyze a certain class of antibody (prob not)
plots = False #if you are running on a local machine and want to see plots - set to True
output_prefix = 'newQuake_post_vax_' #set if you are splitting by barcode and want to prefix the barcode label with something
barcode_file = 'Quake_post_vax_barcodes.txt' #if you are splitting by barcode - input file name of barcode file stored in s3 "info-docs" folder
# format of the barcode file is a text file with one barcode per line:
# barcode1_id <tab> barcode1_sequence
# barcode2_id <tab> barcode2_sequence
# ...etc

min_len = 250 #minimun length of sequences to keep
qual_score = 25 #quality score to keep.  30 is high quality sequences (good for immune rep)
				# 20 or 25 good for most applications


############################################
      ######## Functions ##############
##########################################


def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output


def trim_barcodes(barcode_file, keep_files, output_prefix, data_dir):
	with open(barcode_file) as bf:
		bc_lens={}
		for line in bf:
			bc_lens[line.split('\t')[0]]=len(line.split('\t')[1])

	f_keep_files = []
	for kf in keep_files:
		for key in bc_lens:
			if output_prefix+key == kf:
				print "trimming barcode from "+key
				bash('fastx_trimmer -f %s -i %s -o %s'%(bc_lens[key]+1,data_dir+'/'+kf, data_dir+'/'+kf+'-len-q-t'))
				f_keep_files.append(kf+'-len-q-t')

	return f_keep_files



def clean_seqs(ids, output_prefix='', barcode_file=None, min_len=350, qual_score=30):

	#make dirs to store files
	if not os.path.exists('temp_seqs/'):
		os.makedirs('temp_seqs/')
	else:
		bash('rm temp_seqs/*')

	if not os.path.exists('cleaned_output/'):
		os.makedirs('cleaned_output/')
	else:
		bash('rm cleaned_output/*')


		#copy barcode file
	bash('aws s3 cp s3://info-docs/%s .'%barcode_file)
	outfiles = []
	# loop through files, copy to ec2, process based on input parameters, and copy output back to s3
	for i in ids:
		f = i+'.fastq'

		print "copying %s to ec2"%f
		bash('aws s3 cp s3://patient_repertoire_data/%s .'%f)

		print "removing sequences that are too short"
		bash('fastx_clipper -i %s -o temp_seqs/%s-len.fastq -l %s'%(f,i,min_len))

		print "removing low quality sequences"
		bash('fastq_quality_filter -q %s -p 80 -i temp_seqs/%s-len.fastq -o temp_seqs/%s-len-q.fastq'%(qual_score,i,i))


		if barcode_file:
			print "splitting sequences by barcode"
			split_cmd = 'cat temp_seqs/%s-len-q.fastq | fastx_barcode_splitter.pl --bcfile %s \
						--bol --prefix cleaned_output/%s'%(i,barcode_file,output_prefix)

			keep_files = []
			for line in subprocess.Popen(split_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[1:-3]:
				if int(line.split('\t')[1])>1000:
					keep_files.append(output_prefix+line.split('\t')[0])

			f_keep_files = trim_barcodes(barcode_file, keep_files, output_prefix, 'cleaned_output')


		# copy back to s3
		print "copying %s files back to s3"%i
		for of in f_keep_files:
			bash('fastq_to_fasta -i cleaned_output/%s -o cleaned_output/%s.fasta'%(of,of))
			outfiles.append(of)
			aws_loc = 's3://clean-repertoire-data/' + of +'.fasta'
			bash('aws s3 cp cleaned_output/'+of+'.fasta' + ' '+aws_loc)
		print outfiles


		os.system('rm cleaned_output/*')
		os.system('rm temp_seqs/*')
		os.system('rm %s'%f)

	os.system('rm %s'%barcode_file)
	os.system('rm -r cleaned_output')
	os.system('rm -r temp_seqs')
	return outfiles


def make_dir(dirname):
	if not os.path.exists('%s/'%dirname):
		os.makedirs('%s/'%dirname)
	else:
		os.system('rm %s/*'%dirname)


###########################
####### MAIN SCRIPT #######
###########################

#clean seqs:
# fileIDs = clean_seqs(fileIDs, output_prefix=output_prefix, barcode_file=barcode_file, min_len=min_len, qual_score=qual_score)
fileIDs = ['newQuake_post_vax_MID0-len-q-t', 'newQuake_post_vax_MID1-len-q-t', 'newQuake_post_vax_MID2-len-q-t', \
			'newQuake_post_vax_MID6-len-q-t', 'newQuake_post_vax_MID7-len-q-t','newQuake_post_vax_MID8-len-q-t','newQuake_post_vax_MID9-len-q-t']
# Run VDJ-fasta - store output in s3 folder: ''
multiprocess_vdj_local.run_vdjfasta(fileIDs, num_sequences)


# load in vdj-fasta output and pull out features using rep-seq
ext = '.VDJ.H3.L3.CH1.fa'
s3_dir = 's3://vdjfasta-output/'
counter = 0

for f in fileIDs:
	start = time.time()
	print "started processing", f
	f_string = s3_dir + f + ext

	make_dir('/home/ubuntu/tempvdj')
	make_dir('/home/ubuntu/tree_output')
	make_dir('/home/ubuntu/parsed_fasta')

	#copy from s3 to local
	bash('aws s3 cp %s /home/ubuntu/tempvdj/%s'%(f_string,f+ext))

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



	#remove all temp files that have been transferred
	os.system('sudo rm -f -r /home/ubuntu/tempvdj/*')
	os.system('sudo rm -f -r /home/ubuntu/tree_output/*')
	os.system('sudo rm -f -r /home/ubuntu/parsed_fasta/*')

	counter+=1
	print f, "complete!  %s of %s files finished\n"%(counter,len(fileIDs))
	print "processing time:", time.time()-start

os.system('sudo rm -f -r /home/ubuntu/tempvdj')
os.system('sudo rm -f -r /home/ubuntu/tree_output')
os.system('sudo rm -f -r /home/ubuntu/parsed_fasta')








