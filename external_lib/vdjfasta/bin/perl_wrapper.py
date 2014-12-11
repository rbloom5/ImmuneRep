import subprocess

def runPerl(filedir=''):
	cmd = ["perl", filedir]
	#cmd = ['perl', 'fasta-vdj-pipeline.pl', '--file=/Users/ramantalwar/Documents/Packages/vdjfasta/test/simulated.seqs.fa', '--verbose=1']
	subprocess.call(cmd)

runPerl()