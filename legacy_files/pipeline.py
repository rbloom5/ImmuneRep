import os
import sys
import subprocess
import webbrowser
from Bio import SeqIO

cwd = os.getcwd()

#Set directory paths to SRA toolkit (VDJFasta's paths are set in VDJFasta.pm)

def pipeline(runs):

	#Puts the input runs into list format
	if not isinstance(runs, list):
		runs = [runs]

	for r in runs:

		#Pull in SRA files and convert to .fastq
		print('Converting '+r+' to fastq format...')
		datadir = cwd + '/fastq_files'
		os.system(cwd+'/external_lib/sratoolkit/bin/fastq-dump '+r+' --outdir '+datadir)
		#print 'fastq-dump '+r+' -outdir '+datadir+' '
		print('************************************')

		#Pull in .fastq file and convert to .fasta
		print('Converting '+r+' to fasta formtat...')
		fastq_name = cwd + '/fastq_files/'+r+'.fastq'
		fasta_name = cwd + '/fasta_files/'+r+'.fasta'
		if not os.path.exists(cwd + '/fasta_files'):
			os.makedirs(cwd + '/fasta_files')
		SeqIO.convert(fastq_name, "fastq", fasta_name, "fasta")
		print('************************************')

		#Run FastQC for quality reporting
		qcFileString = cwd + "/QCreports/" + r + "_fastqc.html"

		#Check to see if somebody has already done a QC - if so, just display it
		if os.path.isfile(qcFileString):
			webbrowser.open("file://" + qcFileString)

		#Generate qc report
		else:
			print "Generating QC report..."
			if not os.path.isdir('./QCreports'):
				os.system('mkdir ./QCreports')
			systemString = cwd + "/external_lib/FastQC/fastqc " + fastq_name + " -outdir=./QCreports"
			os.system(systemString)
			webbrowser.open("file://" + qcFileString)
		print('************************************')

		
		#Run VDJFasta
		# print('Running VDJFasta on '+r+'...')
		# cmd = ['perl', cwd + '/external_lib/vdjfasta/bin/fasta-vdj-pipeline.pl', "--file="+fasta_name, "--verbose=1"]
		# subprocess.call(cmd)

pipeline(['SRR1383470'])

