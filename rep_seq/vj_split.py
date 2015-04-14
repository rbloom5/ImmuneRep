#!/usr/bin/python

import re

def write_records(v_j_string,in_file, count, maximum):
    outfile = v_j_string + '.fasta'
    with open(in_file, 'r') as f, open(outfile, 'w') as of:
        for line in f:
            if re.search('^>', line):
                if count >= maximum:
                    return count
                write = 0
                V = line.split(';')[1].split(' ')[0]
                J = line.split(';')[3].split(' ')[0]
                new_v_j_string = V + '-' + J
                if new_v_j_string == v_j_string:
                    write = 1
                    count+=1

            if write:
                of.write(line)

    return count
            
    


def v_j_file_split(filepath, maximum = 10000): #filepath must be to file with XXX.VDJ.H3.L3.CH1.fa from VDJfasta
    with open(filepath, "rU") as in_handle:

        v_j=[]
        count = 0
        for line in in_handle:
            if re.search('^>', line):
                count+=1
                if count > maximum:
                    break
                V = line.split(';')[1].split(' ')[0]
                J = line.split(';')[3].split(' ')[0]
                v_j_string = V + '-' + J
                for i in range(len(v_j_string)):
                    if v_j_string[i] == '/':
                        v_j_string = v_j_string[:i] + '_' + v_j_string[i+1:]

                if v_j_string not in v_j:
                    v_j.append(v_j_string)
                    count = write_records(v_j_string, filepath, count, maximum)
            


def v_split(Reads):
    #Rep must be dict of Ab_read objects
    Vs = []
    Vdict = {}
    searchable_Reads = Reads.copy() 

    for item in Reads:
        if Reads[item].V:
            currentV = Reads[item].V[0] #if there is more than 1 V match, arbitrarily picks first one
            if currentV not in Vs:
                Vs.append(currentV)
                Vdict[currentV] = {}

                IDs = [] #this is a dict where I store the keys that match currentV so
                        #I can delete them afterwards and not waste search time
                        #iterating over them again
                for item2 in searchable_Reads:
                    if searchable_Reads[item2].V:
                        if searchable_Reads[item2].V[0] == currentV:
                            IDs.append(item2)
                            Vdict[currentV][item2]  = searchable_Reads[item2]
                for ID in IDs:
                    del searchable_Reads[ID]

    return Vdict


def v_first4_split(Reads):
    #Rep must be dict of Ab_read objects
    #this groups Rep into V segments based on first four characters (i.e. IVH4)
    # rather than the whole thing (i.e. IVH4-33)
    Vs = []
    Vdict = {}
    searchable_Reads = Reads.copy() 

    for item in Reads:
        if Reads[item].V:
            currentV = Reads[item].V[0][0:5] #if there is more than 1 V match, arbitrarily picks first one
            if currentV not in Vs:
                Vs.append(currentV)
                Vdict[currentV] = {}

                IDs = [] #this is a dict where I store the keys that match currentV so
                        #I can delete them afterwards and not waste search time
                        #iterating over them again
                for item2 in searchable_Reads:
                    if searchable_Reads[item2].V:
                        if searchable_Reads[item2].V[0][0:5] == currentV:
                            IDs.append(item2)
                            Vdict[currentV][item2]  = searchable_Reads[item2]
                for ID in IDs:
                    del searchable_Reads[ID]

    return Vdict








#v_j_file_split('SRR1298742.VDJ.H3.L3.CH1.fa')


