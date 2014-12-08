import sys
import subprocess

def run_VDJFasta():
	fileName = "--file=" + sys.argv[1]
	progName = sys.argv[2]
	#cmd = ["perl","fasta-vdj-pipeline.pl",fileName,"--verbose=1"]
	cmd = ["perl",progName,fileName,"--verbose=1"]
	subprocess.call(cmd)

run_VDJFasta()