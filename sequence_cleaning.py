#!/usr/bin/python

import os
import subprocess

def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output

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


	# loop through files, copy to ec2, process based on input parameters, and copy output back to s3
	for i in ids:
		f = i+'.fastq'

		print "copying %s to ec2"%f
		bash('aws cp s3://patient_repertoire_data/%s .'%f)

		print "removing sequences that are too short"
		bash('fastx_clipper -i %s -o temp_seqs/%s-len.fastq -l %s',%(f,i,min_len))

		print "removing low quality sequences"
		bash('fastq_quality_filter -q %s -p 80 -i temp_seqs/%s-len.fastq -o temp_seqs/%s-len-q.fastq',(qual_score,i,i))


		if barcode_file:
			print "splitting sequences by barcode"
			bash('cat temp_seqs/%s-len-q.fastq | fastx_barcode_splitter.pl \
				--bcfile %s --bol --exact --prefix cleaned_output/%s',(i,barcode_file,output_prefix))

		# copy back to s3
		print "copying %s files back to s3"%i
		output_files = [i for i in os.listdir('cleaned_output/')]
		for of in output_files:
			aws_loc = 's3://clean-repertoire-data/' + of
			bash('aws s3 cp cleaned_output/'+of + ' '+aws_loc)


		os.system('rm cleaned_output/*')
		os.system('rm temp_seqs/*')
		os.system('rm %s',f)









