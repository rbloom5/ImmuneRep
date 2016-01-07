#!/usr/bin/python

import os
import subprocess
import gzip

def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output

def make_dir(dirname):
	if not os.path.exists('%s/'%dirname):
		os.makedirs('%s/'%dirname)
	else:
		os.system('rm %s/*'%dirname)

def copy_to_s3(files, local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://%s/'%s3_dir + f
		bash('aws s3 cp %s/'%local_dir+f+' '+aws_loc)

def copy_from_s3(files,local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://%s/'%s3_dir + f
		bash('aws s3 cp ' + aws_loc + ' %s/'%local_dir + f)



def d50(metadata, diversity_file, label_column=None, group_column=None):
	#that specifies which group the sample goes into (i.e. control vs disease)

	infiles = []
	labels=[]
	groups=[]
	# ages=[]
	with open(metadata) as f:
		f.readline()
		for line in f:
			infiles.append(line.split('\t')[0])
			if label_column:
				labels.append(line.strip().split('\t')[label_column])
			if group_column:
				groups.append(line.strip().split('\t')[group_column])
			# ages.append(line.split('\t')[2])

	with open(diversity_file, 'w') as of:

		of.write('sample'+'\t' + 'd50' + '\t' + 'label' + '\t' + 'group' + '\n' )
		for i in range(len(infiles)):
			infile = infiles[i]
			# os.system('gunzip -k %s'%infile)

			with open(infile) as f:
				d50_not_reached=1
				d50_clone=0
				clone_count=0
				read_count=0
				total_clones=0
				f.readline()
				for line in f:
					total_clones+=1
					read_count+=float(line.strip().split('\t')[1])
					clone_count+=1
					if read_count>=.5 and d50_not_reached:
						d50_clone=clone_count
						d50_not_reached=0
			# os.system('rm %s'%infile)
			if label_column and group_column:
				of.write(infile + '\t' + str(d50_clone/float(total_clones)) + '\t' + labels[i] + '\t' + groups[i] + '\n')
			else:
				of.write(infile + '\t' + str(d50_clone/float(total_clones)) + '\n')
			print "done with ",infile




def fastq_dump(f, clean=True):
	try:
		print "downloading %s from SRA"%ids
		os.system('fastq-dump %s'%f)
		copy_to_s3(f, '.', 't-cell-rep/raw-tcr-seqs')

	except:
		return 'error getting %s in fastq-dump'%ids



def run_format_MITCR(ids, meta_file = 'metadata.txt'):

	# print "running MiTCR on %s"%ids
	try:
		print 'java -jar ./software/mitcr.jar -pset jprimer %s %s-mitcr.txt'%(ids, ids[:-6])
		print bash('java -jar ./software/mitcr.jar -pset jprimer %s %s-mitcr.txt'%(ids, ids[:-6]))
		print 'java -jar ./software/vdjtools-1.0-SNAPSHOT.jar Convert %s-mitcr.txt formatted'%ids[:-6]
		print bash('java -jar ./software/vdjtools-1.0-SNAPSHOT.jar Convert %s-mitcr.txt formatted'%ids)
		copy_to_s3('formatted.%s-mitcr.txt'%ids.split('/')[-1][:-6], '.', 't-cell-rep/mitcr-output')
		# os.system('rm %s.fastq'%ids)
		os.system('rm %s'%ids)
		with open(meta_file,'a') as f:
			f.write('formatted.%s\n'%ids)
			os.system('rm formatted.%s-mitcr.txt'%ids)

	except:
		print 'error running MITCR on %s'%ids



def get_stats(metadata, output_prefix='', plots=True, factor='', label=''):
	print 'getting rep statistics'

	#basic statistics
	status1 = bash('java -jar ./software/vdjtools-1.0-SNAPSHOT.jar CalcBasicStats -m %s /home/ubuntu/MITCR_out/%s '%(metadata,output_prefix))
	print status1
	if status1:
		print "error getting basic statistics"



	# get and plot V and J usage
	seg_string = 'java -jar ./software/vdjtools-1.0-SNAPSHOT.jar CalcSegmentUsage '
	print seg_string
	if plots:
		seg_string+='-p '
	if factor:
		seg_string+='-f %s '%factor
	if label:
		seg_string+='-l %s '%label

	seg_string+= '-m %s /home/ubuntu/MITCR_out/%s '%(metadata,output_prefix)
	print seg_string
	status2=bash(seg_string)
	print status2
	if status2:
		print "error getting segment usage"


	#Get Spectratype (length of cdr3 region)
	print 'java -jar ./software/vdjtools-1.0-SNAPSHOT.jar CalcSpectratype -m %s /home/ubuntu/MITCR_out/%s '%(metadata,output_prefix)
	status3=bash('java -jar ./software/vdjtools-1.0-SNAPSHOT.jar CalcSpectratype -m %s /home/ubuntu/MITCR_out/%s '%(metadata,output_prefix))
	print status3
	if status3:
		print "error getting spectractype"	


	##Get diversity stats
	print 'java -jar ./software/vdjtools-1.0-SNAPSHOT.jar CalcDiversityStats -m %s /home/ubuntu/MITCR_out/%s '%(metadata,output_prefix)
	status4=bash('java -jar ./software/vdjtools-1.0-SNAPSHOT.jar CalcDiversityStats -m %s /home/ubuntu/MITCR_out/%s '%(metadata,output_prefix))
	print status4
	if status4:
		print "error getting diversity"


	## Get Rarefaction plot
	if plots:
		rare_string = 'java -jar ./software/vdjtools-1.0-SNAPSHOT.jar RarefactionPlot '
		if factor:
			rare_string+='-f %s '%factor
		if label:
			rare_string+='-l %s '%label
	rare_string+= '-m %s /home/ubuntu/MITCR_out/%s '%(metadata,output_prefix)
	print rare_string			
	status5=os.system(rare_string)
	if status5:
		print "error calculating rarefaction plot"


	#get quantile plots
	infiles=[]
	if plots:
		with open(metadata) as f:
			f.readline()
			for line in f:
				infiles.append(line.split('\t')[0])
		for inf in infiles:
			head, tail = os.path.split(inf)
			status6=bash('java -jar ./software/vdjtools-1.0-SNAPSHOT.jar PlotQuantileStats %s /home/ubuntu/MITCR_out/quantile_plots/%s '%(inf,tail))
			if status6:
				print "error plotting quantile plots"





# TCR_files = [\

			
# 			]

# MS_files = ['SRR1956786-mitcr.txt',	\
# 'SRR1956783-mitcr.txt',\
# 'SRR1956784-mitcr.txt',	\
# 'SRR1956785-mitcr.txt',	\
# 'SRR1956788-mitcr.txt',	\
# 'SRR1956790-mitcr.txt',	\
# 'SRR1956789-mitcr.txt',\
# 'SRR1956787-mitcr.txt',\
# 'SRR1956792-mitcr.txt',\
# 'SRR1956791-mitcr.txt',\
# 'SRR1956794-mitcr.txt',\
# 'SRR1956793-mitcr.txt',\
# 'SRR1956796-mitcr.txt',\
# 'SRR1956795-mitcr.txt',]

ITP_files = ['TPO-113.assembled',\
			'TPO-114.assembled',\
			'TPO-118.assembled',\
			'TPO-120.assembled',\
			'TPO-137.assembled',\
			'TPO-139.assembled',\
			'TPO-92.assembled',\
			'TPO-95.assembled',\
			]

import glob

# TCR_files = glob.glob('/home/ubuntu/sequence_files/Zehnder_RNA/*.fastq')
# for ids in ITP_files:
# 	i=ids+'.fastq'
# 	copy_from_s3(i, '.', 't-cell-rep/mitcr-output')
# 	run_format_MITCR(i, meta_file = 'TPO_ITP_meta.txt')

get_stats('zehnder_metadata.txt', output_prefix='zehnder_all', plots=True, factor='disease', label='sample.id')
d50('zehnder_metadata.txt','MITCR_out/zehnder_all_d50.txt', label_column=1, group_column=2)































