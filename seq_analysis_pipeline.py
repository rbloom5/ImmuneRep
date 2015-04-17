#!/usr/bin/env python 

import Rep_sequence_analysis
reload(Rep_sequence_analysis)

dir = '/Users/Ryan/SoftwareProjects/Databases/fastq/'

files = [dir+'SRR1298742.VDJ.H3.L3.CH1.fa',\
            dir+'SRR1383326.VDJ.H3.L3.CH1.fa',\
            dir+'SRR1383453.VDJ.H3.L3.CH1.fa',\
            dir+'SRR1383463.VDJ.H3.L3.CH1.fa',\
            dir+'SRR1383472.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_0.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_1.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_2.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_3.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_6.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_7.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_8.VDJ.H3.L3.CH1.fa',\
#             dir+'QuakeMID_9.VDJ.H3.L3.CH1.fa',\
            ]

RepsSmall = Rep_sequence_analysis.Rep_seq(files,num=10000)
print "loaded"
RepsSmall.find_clones()
print "finding clones done"
RepsSmall.find_clusters()

