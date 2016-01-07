 #!/usr/bin/python

import gzip
import os

infiles = []
ids=[]
ages=[]
with open('all_C_metadata.txt') as f:
    f.readline()
    f.readline()
    for line in f:
        infiles.append(line.split('\t')[0])
        ids.append(line.split('\t')[1])
        ages.append(line.split('\t')[2])

with open('all_C_samples/diversity.txt', 'w') as of:
    of.write('sample'+'\t' + 'age' + '\t' + 'd50' + '\n')
    for i in range(len(infiles)):
        infile = infiles[i]
        os.system('gunzip -k %s'%infile)

        with open(infile[:-3]) as f:
            print infile
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
        os.system('rm %s'%infile[:-3])
        of.write(ids[i] + '\t' + ages[i] + '\t' + str(d50_clone/float(total_clones))+'\n')