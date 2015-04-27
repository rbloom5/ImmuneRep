import os
from Bio import SeqIO
import subprocess
import re
import json
import time

#filenames=[
#'SRR1383450',
#'SRR1383455',
#'SRR1383466',
#'SRR1383470',]
# 'SRR1383448',
# 'SRR1383472',
# 'SRR1383473',
# 'SRR1383474']
#'SRR747232',
#'SRR747758',
#'SRR747760',
#'SRR747785',
#'SRR747766',
#'SRR747768',
#'SRR765688']
# filenames = ['SRR1383463']



def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output

def create_chunks(filenames, data_dir):
    #filenames should not have .fasta extension - just SRR+numbers
    #formatting inputs
    if data_dir[-1]!='/':
        data_dir = data_dir + '/'
    
    if not isinstance(filenames, list):
        filenames = [filenames]
        
        
        
    #list of all the filenames created    
    all_files = []
    
    #loop through filenames    
    for f in filenames:
        #formatting - adding .fasta extension to the input strings

        file_string = data_dir + f + '.fasta'

    
    
    
        #open up the fasta file
        with open(file_string, "rU") as in_handle:
            records = 0
            out_files = 1
            out_file_string = data_dir + f + '_' + str(out_files)+'.fasta'
            all_files.append(f + '_' + str(out_files))
            out_handle = open(out_file_string, 'w+')
            
            #loop through lines in file and count the number of '>' (beginning of each fasta record)
            for line in in_handle:
                if re.search('^>', line):
                    records+=1
                
                #if you get to 2500 records, close the old file, open up a new file and start writing to that
                if records == 1000:
                    records = 0
                    out_files +=1
                    out_handle.close()
                    out_file_string = data_dir + f + '_' + str(out_files)+'.fasta'
                    all_files.append(f + '_' + str(out_files))
                    out_handle = open(out_file_string, 'w+')
                
                
                #write the line to the out file    
                out_handle.write(line)
            out_handle.close()
        os.remove(file_string)
    return all_files
            
    
#create_chunks(['SRR1383460'], '/Users/Ryan/SoftwareProjects/ImmuneRep/fasta_files')

def merge_chunks(chunknames,filename, datadir,extensions):
#filenames should simply be 'SRR' followed by numbers - no extension

    #format inputs
    if datadir[-1]!='/':
        datadir = datadir + '/'
    
    for ext in extensions:
	    main_string = datadir+filename+ext
	    with open(main_string, 'w+') as main_file:
		    for chunkname in chunknames:
			    with open(datadir+chunkname+ext,'r') as chunk_file:
				    for line in chunk_file:
					    if line[:13]!='WRITEDONEFLAG':
						    main_file.write(line)
			    bash('rm -f '+datadir+chunkname+ext)

def humanize_time(secs):
    mins, secs = divmod(secs, 60)
    hours, mins = divmod(mins, 60)
    return '%02d:%02d:%02d' % (hours, mins, secs)


def run_vdjfasta(filenames):

	extensions = ('.C.germdata.txt','.L1.acc.txt','.aa.fa','.wIgs.stock',\
              '.D.germdata.txt','.L2.acc.txt',\
              '.wIgs.unique.fa','.H1.acc.txt','.L3.acc.txt','.fasta',\
              '.H2.acc.txt','.V.germdata.txt','.wIgs.Vh-gs-Vk.a2m',\
              '.J.germdata.txt','.H3.acc.txt','.VDJ.H3.L3.CH1.fa',\
              '.wIgs.Vh-gs-Vk.c2m','.VDJ.coords.txt','.wIgs.fa','.aa.fa.1e-10.score.txt')

	s3dir='patient_repertoire_data/'
	datadir='/home/ubuntu/data/'
	Nmax=32

	runbase='sudo docker run -d -v /home/ubuntu/data:/home/data rzeller/vdj-process /home/vdjfasta/bin/fasta-vdj-pipeline.pl --file=/home/data/'

	existre=re.compile(r'Error response from daemon: Cannot start container (.*): file exists')
	processtime=time.time()
	for filename in filenames:
		try:
			filetime=time.time()
			#copy .fastq from s3
			# cmd = 'aws s3 cp s3://'+s3dir+filename+'.fastq '+datadir+filename+'.fastq'
			# os.system(cmd)

			# # convert to .fasta file
			# SeqIO.convert(datadir+filename+'.fastq', 'fastq', datadir+filename+'.fasta', 'fasta')

			# # copy .fasta file to s3
			# cmd='aws s3 cp '+datadir+filename+'.fasta s3://'+s3dir+filename+'.fasta'
			# os.system(cmd)

			# # delete .fastq from ec2
			# cmd = 'rm '+datadir+filename+'.fastq'
			# os.system(cmd)

			cmd = 'aws s3 cp s3://'+s3dir+filename+'.fasta '+datadir+filename+'.fasta'
			os.system(cmd)

			# break into chunks
			chunknames=create_chunks(filename,datadir)
			Nchunks=len(chunknames)

			# initialize chunks status dictionary
			chunks={}
			for chunkname in chunknames:
				chunks[chunkname]={'status':2}

			flag=1
			containers=[]
			while flag:
				Nrunning=0
				runningchunks=[key for key in chunks.keys() if chunks[key]['status']==1]
				for chunkname in runningchunks:
					running=json.loads(bash('sudo docker inspect '+chunks[chunkname]['container']))[0]['State']['Running']
					if running:
						Nrunning+=1
					else:
						if len(bash('cat '+datadir+chunkname+'.VDJ.H3.L3.CH1.fa').split('\n'))>10:
							chunks[chunkname]['status']=0
							bash('sudo docker rm '+chunks[chunkname]['container'])
							# bash('''sudo docker ps -a | grep Exit | awk '{print $1}' | xargs --no-run-if-empty sudo docker rm''')
							Ndone=len([key for key in chunks.keys() if chunks[key]['status']==0])
							print repr(Ndone)+'/'+repr(Nchunks),
							print 'estimated time remaining '+humanize_time((time.time()-filetime)/Ndone*(Nchunks-Ndone))
							# print '...finished chunk '+chunkname
						elif (time.time()-chunks[chunkname]['starttime'])>60:
							print 'processing failed on '+filename+'!'
							print 'container '+chunks[chunkname]['container']+' failed while processing '+chunkname
							print bash('sudo docker logs '+chunks[chunkname]['container'])
							bash('sudo docker rm '+chunks[chunkname]['container'])
							bash('''sudo docker ps | grep 'ago' | awk '{print $1}' | xargs --no-run-if-empty sudo docker stop''')
							bash('''sudo docker ps -a | grep 'ago' | awk '{print $1}' | xargs --no-run-if-empty sudo docker rm''')
							break
						else:
							pass

				time.sleep(.5)
				Nopen=Nmax-Nrunning
				rawchunks=[key for key in chunks.keys() if chunks[key]['status']==2]
				Ncreate=min(Nopen,len(rawchunks))
				for chunkname in rawchunks[:Ncreate]:
					try:
						output=bash(runbase+chunkname+'.fasta'+' --verbose=1').split('\n')[0]
						chunks[chunkname]['container']=output
						chunks[chunkname]['status']=1
						chunks[chunkname]['starttime']=time.time()
					except:
						pass
						# bash('''sudo docker ps -a | grep Exit | awk '{print $1}' | xargs --no-run-if-empty sudo docker rm''')

				if not Nrunning+len(rawchunks):
					flag=0
				else:
					# bash('''sudo docker ps -a | grep Exit | awk '{print $1}' | xargs --no-run-if-empty sudo docker rm''')
					time.sleep(.5)

			merge_chunks(chunknames,filename,datadir,extensions)

			for ext in extensions:
				bash('aws s3 cp '+datadir+filename+ext+' s3://vdjfasta-output/')
				bash('rm -f '+datadir+filename+ext)

			print 'processing of ' +filename+ 'completed!'
			print 'file completion time'+humanize_time(time.time()-filetime)
		except:
			bash('''sudo docker ps | grep 'ago' | awk '{print $1}' | xargs --no-run-if-empty sudo docker stop''')
			bash('''sudo docker ps -a | grep 'ago' | awk '{print $1}' | xargs --no-run-if-empty sudo docker rm''')
			bash('rm -f data/*')				

	print 'Processing Complete!'
	print 'process completion time'+humanize_time(time.time()-processtime)



		

		 





