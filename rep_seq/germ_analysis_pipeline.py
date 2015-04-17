#!/usr/bin/python

import Rep_sequence_analysis

dir = '/Users/Ryan/SoftwareProjects/ImmuneRep/fastq_files/'

MSfiles = [dir+'SRR1298742.VDJ.H3.L3.CH1.fa',\
			dir+'SRR1383326.VDJ.H3.L3.CH1.fa',\
			dir+'SRR1383453.VDJ.H3.L3.CH1.fa',\
			dir+'SRR1383463.VDJ.H3.L3.CH1.fa',\
			dir+'SRR1383472.VDJ.H3.L3.CH1.fa',\]


Controlfiles = [dir+'QuakeMID_0.VDJ.H3.L3.CH1.fa',\
				dir+'QuakeMID_1.VDJ.H3.L3.CH1.fa',\
				dir+'QuakeMID_2.VDJ.H3.L3.CH1.fa',\
				dir+'QuakeMID_3.VDJ.H3.L3.CH1.fa',\
				dir+'QuakeMID_6.VDJ.H3.L3.CH1.fa',\
				dir+'QuakeMID_7.VDJ.H3.L3.CH1.fa',\
				dir+'QuakeMID_8.VDJ.H3.L3.CH1.fa',\
				dir+'QuakeMID_9.VDJ.H3.L3.CH1.fa',\]

MS_V_freqs = {}
MS_f4V_freqs = {}
MS_top_clone = {}
MS_top_ten_clones = {}

for f in MSfiles:
	currentRep = Rep_sequence_analysis.Rep_seq(f, 50000)
	currentRep.split_V_4()
	currentRep.find_clones()

	MS_V_freqs[f] = currentRep.V_fractions
	MS_f4V_freqs[f] = currentRep.f4V_fractions
	MS_top_clone[f] = currentRep.Clones[1].percent_reads

	sum10 = 0
	for i in range(1,11):
		sum10 += currentRep.All_clones[i].percent_reads
	MS_top_ten_clones[f] = sum10
	print f + " complete!"



C_V_freqs = {}
C_f4V_freqs = {}
C_top_clone = {}
C_top_ten_clones = {}

for f in Controlfiles:
	currentRep = Rep_sequence_analysis.Rep_seq(f, 50000)
	currentRep.split_V_4()
	currentRep.find_clones()

	C_V_freqs[f] = currentRep.V_fractions
	C_f4V_freqs[f] = currentRep.f4V_fractions
	C_top_clone[f] = currentRep.Clones[1].percent_reads

	sum10 = 0
	for i in range(1,11):
		sum10 += currentRep.All_clones[i].percent_reads
	C_top_ten_clones[f] = sum10
	print f + " complete!"






