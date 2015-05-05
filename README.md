ImmuneRep
=========
For instructions on loading data into python, see the data retrieval page in the wiki.

packages required for rep-seq:
multiprocessing, 
joblib,
sklearn,
numpy,
Biopython,
seaborn,
ete,
aws command line tools,
pandas,


all of which can be obtained through pip/easy-install etc.  

the classes used within the Rep-seq class (AB_read, Clone, Cluster) are in the file rep-classes.py.  

initialize:
Rep_seq(filepath, num=False, ABtype=[])

the rep seq class is made to hold and calculate all the useful properties from a patients b-cell repertoire

inputs:

file path: the path to the .VDJ.H3.L3.CH1.fa from vjd-fasta
num: the maximum number of reads you want to analyze.  If false, it will take all the reads in the file
ABtype: if you only want to take antibodies of a certain type (i.e. IGHG, IGHM etc.)  if blank, it will take all the antibodies regardless of type

attributes created:
.Reads: a dict of AB_read objects, the key to the dict is the sequence id number from vdj-fasta output
.Reads_split_by_VJ: reads organized by their V and J gremlins used.  format: Reads_split_by_VJ = {VJgermline :{ID: AB_read object}}
.num_Reads: total number of reads


.split_V_4():
method that splits reads by the major V segment (IGV1, IGV2 etc.) and calculates frequencies.

returns: nothing
attributes created: 
.Reads_split_by_f4V: a dict of the form Reads_split_by_f4V = {IGVH1:{ID: Ab read object}…}
.VJ_freqs: a dict that tells you how many reads have the corresponding V segment - form:  VJ_freqs={IGVH1: number of reads, …}
.VJ_fractions: same as frees, but instead of number of reads it is the fraction of total reads



.find_clones(parallel=False)
This creates a dict of clone objects from all the sequences in reads a clone object combines sequences that are identical as well as within a Hamming distance of 1 from each other (to account for sequencing errors - see von Budingen, 2012)

if you want it to use multiple cores in parallel to perform the clustering, set parallel=True

returns: nothing
attributes created:
.Clones: a dict where the keys are number that correspond to the read rank of the clone (i.e. the close that appears the most often in the repertoire has the key ‘1’)
.Clones_split_by_VJ: clones organized by the VJ gremlin segments they use.  form is a dict of dicts, same form as Reads_split_by_VJ
.num_clones: number of unique clones



.find_clusters()


.get_stats():
calculates useful statistics from the repertoire

returns: nothing
attributes created:
.sh_dict: a dict of somatic hypermutation distributions.  The keys of the dict are the ab type (IGHG, IGHM etc.) the value is a list of the number of somatic hypermutations for each clone of that type.  
.median_sh: a dict with ab type as the keys, where the value is the median amount of somatic hypermutation from clones of that ab type
.clone_distribution: a list 100 items long.  each value is the percent of total reads that the percentile of clones takes up.  (i.e. if the top 1% of clones make up 20% of the total reads in the repertoire, then the first value in the list will be .2).  Useful for seeing if particular clones are over represented compared to the norm
.top_clone_fraction: the fraction of the repertoire that just the top clone takes up
.top_10_fraction: the fraction of the repertoire that the top 10 clones take up
.top_100_fraction: the fraction of the repertoire that the top 100 clones take up

.ABtype_fractions: a dict where the keys are the ab type and the value is the fraction of the total reads in the repertoire that is made up of antibodies of this type
.ABtype_unique: a dict where keys are ab type and the value is the fraction of unique reads (instead of total reads) that are made up of antibodies of this type



.plots()

makes a bunch of cool plots.  If you have seaborn, the plots will be prettier


