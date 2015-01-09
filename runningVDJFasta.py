import sys
import subprocess
from Bio import SeqIO

#This script acts as a python wrapper to the perl program VDJFasta. 
#VDJFasta is a script that takes in a .fastq file of VDJ sequenced regions and outputs a 
#.fasta containting the V, J regions and the amino acid sequence for the CDR3

#This program requires you to download and setup VDJFasta. Also you need to point to your location of 
#fasta-vdj-pipeline.pl when you call the prgoram as the second argument

#The format to run runningVDJFasta.py is...
#runningVDJFasta.py /location/sequence.data.fastq /location/fasta-vdj-fasta-pipeline.pl


def run_VDJFasta():
	#This is for testing purposes
	#fastqName = "/Users/ramantalwar/Documents/Python/fasta_conversion_testing/SRR1298742.fastq"
	#progName = "/Users/ramantalwar/Documents/Packages/vdjfasta/bin/fasta-vdj-pipeline.pl"

	#pulls in the file and coverts it from .fastq to .fasta
	fastqName = sys.argv[1]
	fastaName = "testing.fasta"
	SeqIO.convert(fastqName, "fastq", fastaName, "fasta")
	fileName = "--file=" + fastaName

	#pulls in location of vdj-fasta-pipeline.pl
	progName = sys.argv[2]

	cmd = ["perl",progName,fileName,"--verbose=1"]
	subprocess.call(cmd)

run_VDJFasta()
