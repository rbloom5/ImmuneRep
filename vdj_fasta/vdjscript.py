#!/usr/bin/python
import subprocess
from subprocess import PIPE
import os

def cmdline(command):
    process = subprocess.Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]

cwd = os.getcwd()
# for i in range(1,6):
# 	fasta_name  = cwd+'/fastq_files/QuakeMID_3_' + str(i) + ".fasta"
# 	cmd = 'perl ' + cwd + "/external_lib/vdjfasta/bin/fasta-vdj-pipeline.pl --file="+fasta_name, "--verbose=1"
# 	out = cmdline(cmd)
# 	print out

# for i in range(1, 10):
# 	fasta_name  = cwd+'/fastq_files/QuakeMID_6_' + str(i)+ ".fasta"
# 	cmd = 'perl ' + cwd + "/external_lib/vdjfasta/bin/fasta-vdj-pipeline.pl --file="+fasta_name, "--verbose=1"
# 	out = cmdline(cmd)
# 	print out

# for i in range(1, 7):
# 	fasta_name  = cwd+'/fastq_files/QuakeMID_7_' + str(i)+ ".fasta"
# 	cmd = 'perl ' + cwd + "/external_lib/vdjfasta/bin/fasta-vdj-pipeline.pl --file="+fasta_name, "--verbose=1"	
# 	print cmd
# 	out = cmdline(cmd)
# 	print out

for i in range(1, 4):
	fasta_name  = cwd+'/fastq_files/QuakeMID_8_' + str(i)+ ".fasta"
	cmd = 'perl ' + cwd + "/external_lib/vdjfasta/bin/fasta-vdj-pipeline.pl --file="+fasta_name, "--verbose=1"
	out = cmdline(cmd)
	print out

for i in range(1, 7):
	fasta_name  = cwd+'/fastq_files/QuakeMID_9_' + str(i)+ ".fasta"
	cmd = 'perl ' + cwd + "/external_lib/vdjfasta/bin/fasta-vdj-pipeline.pl --file="+fasta_name, "--verbose=1"	
	out = cmdline(cmd)
	print out



