import os
import sys
import subprocess
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
		SeqIO.convert(fastq_name, "fastq", fasta_name, "fasta")
		print('************************************')

		#Run VDJFasta
		print('Running VDJFasta on '+r+'...')
		cmd = ['perl', cwd + '/external_lib/vdjfasta/bin/fasta-vdj-pipeline.pl', "--file="+fasta_name, "--verbose=1"]
		subprocess.call(cmd)

pipeline('SRR735691')

