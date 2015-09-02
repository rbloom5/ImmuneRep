#!/usr/bin/python

import os
import glob


complete = len(os.listdir('/home/ubuntu/parsed_fasta'))/2
file_dir = [i for i in glob.glob('/home/ubuntu') if i.endswith('_vj_files')]
f= glob.glob('/home/ubuntu/*_vj_files')
total = len(os.listdir(f[0]))

print '%s out of %s VJ pairs complete'%(complete,total)