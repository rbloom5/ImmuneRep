#! /usr/bin/env python

import sys
import os
import webbrowser
import subprocess


########  DataQC.py ########
# this file contains 3 functions:
# QCreport: generates and displays the sequencing data quality report using FastQC
# getdata: downloads and converts .sra files into .fastq format
# 


def QCreport(filename):
	# generates and displays sequencing quality reports from FastQC
	# filename: the path of the .fastq file that you want to analyze

	qcFileString = os.path.dirname(os.path.realpath(__file__)) + "/QCreports/" + filename[-16:-6] + "_fastqc.html"

	#Check to see if somebody has already done a QC - if so, just display it
	if os.path.isfile(qcFileString):
		webbrowser.open("file://" + qcFileString)

	#Generate qc report
	else:
		print "Generating QC report..."
		systemString = os.path.dirname(os.path.realpath(__file__)) + "/FastQC/fastqc " + filename + " -outdir=./QCreports"
		os.system(systemString)
		webbrowser.open("file://" + qcFileString)


#QCreport("/Users/Ryan/PythonProjects/Databases/ImmRep/Datasets/SRR1298742.fastq")


def getdata(runs, datadir = ''):
	# downloads and converts .sra files to fastq format

	#### inputs #####
	# runs: the run numbers (as a string or list) for the .sra files you want (i.e. 'SRR...').
	#       The script will automatically find, download, and convert the files

	# datadir = the directory where you want the fastq files to go

	#make list so we can iterate through
	if not isinstance(runs, list):
		runs = [runs]

	for r in runs:
		print('converting '+r+' to fastq format...')
		os.system('fastq-dump '+r+ ' --outdir '+datadir)
		print 'fastq-dump '+ r + ' -outdir '+datadir+' '
		print('************************************')


			#### OLD ####

			# if len(sys.argv)<2:
			# 	datadir=''
			# else:
			# 	datadir=sys.argv[1]


			#example of full ftp path:
			#ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR138/SRR1383465

			# ftpbase='ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'
					# ftpstring=ftpbase+r[:6]+'/'+r+'/'+r+'.sra'
					# filestring=datadir+r+'.sra'
					# fastqfilestring=datadir+r+'.fastq'

					# #download file
					# print('downloading '+r+'...')
					# urllib.urlretrieve(ftpstring,filestring)

					#convert to .fastq



#getdata(runs, datadir='~/PythonProjects/Databases/ImmRep/Datasets')


def getdata_QC(runs, datadir=''):
	# Combines getdata and dataQC
	# will download and converts .sra files to fastq format then generate and display QC

	#### inputs #####
	# runs: the run numbers (as a string or list) for the .sra files you want (i.e. 'SRR...').
	#       The script will automatically find, download, and convert the files

	# datadir = the directory where you want the fastq files to go, default is directory of this script


	if datadir == '':
		datadir = os.path.dirname(os.path.realpath(__file__))

	#make runs a list if it is not already
	if not isinstance(runs, list):
		runs = [runs]

	for r in runs:
		#download and convert data to fastq format
		getdata(r, datadir)

		#This just adds some flexibility depending on whether the user 
		#   put a '/' on the end of the dirname or not
		if datadir[-1] == '/':
			fastqfilestring=datadir+r+'.fastq'
		else:
			fastqfilestring=datadir+'/' + r+'.fastq'

		#generate qc report
		print "Generating QC report..."
		systemString = os.path.dirname(os.path.realpath(__file__)) + "/FastQC/fastqc " + fastqfilestring + " -outdir=./QCreports"
		os.system(systemString)

		qcFileString = os.path.dirname(os.path.realpath(__file__)) + "/QCreports/" + r + "_fastqc.html"
		webbrowser.open("file://" + qcFileString)


runs=['SRR1298742', 'SRR1298742']
# 'SRR1383465','SRR1383464','SRR1383463',
# 'SRR1383462','SRR1383461','SRR1383460',
# 'SRR1383453','SRR1383452','SRR1383451',
# 'SRR1383448','SRR1383447','SRR1383326',
# 'SRR1383474','SRR1383473','SRR1383472']

# getdata_QC(runs, '~/PythonProjects/Databases/ImmRep/Datasets')
