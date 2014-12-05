#! /usr/bin/env python

import sys
import os
import webbrowser

def DataVis(filename):
	systemString = "FastQC/fastqc " + filename + " -outdir=./QCreports"
	os.system(systemString)
	cwd = os.getcwd()
	htmlString = "file://" + cwd + "/QCreports/" + filename[-16:-6] + "_fastqc.html"
	webbrowser.open(htmlString)


#DataVis("/Users/Ryan/PythonProjects/Databases/ImmRep/Datasets/SRR1383461.fastq")
