
This software performs analysis of antibody variable domains. The principles behind the software were described in 
# Glanville J, Zhai W, Berka J et al. Precise determination of the diversity 
# of a combinatorial antibody library gives insight into the human immunoglobulin 
# repertoire. Proc Natl Acad Sci USA. 2009;106:20216–20221
  
The software was released as part of the methods of that publication, and is made available to the community 
for non-commercial purposes. The software is provided "as is," with no implied support. For those interested in
contributing code to the package, the authors can be reached at Jacob <dot> Glanville <at> pfizer <dot> com.

Dependencies
============
VDJFasta requires NCBI's blast toolkit, HMMER's hmm toolkit, and Perl. It is known to work correctly with

* ncbi-blast-2.2.25+
* hmmer-3.0 (Version 3.0; March 2010)
 
Installation
============

The package can be unpacked as
 
$ tar -xzvf vdjfasta-1.0.tgz

This will create the following directory structure

vdjfasta/bin/			<-- bin directory with executable perl scripts
        /db/			<-- blast and hmm reference databases
        /lib/			<-- perl module libraries
        /test/			<-- example sequences to test the software
        /README.txt		<-- this readme file

Once unpacked, a number of paths will need to be set in the constructor of vdjfasta/lib/VDJFasta.pm
Adjust the paths to your current Blast and HMMER installations:

  $self->{blast}        = "/tools/apps/NCBI/blast-2.2.17/bin/blastall";
  $self->{hmmsearch}    = "/tools/apps/HMMER/hmmsearch";
  $self->{hmmalign}     = "/tools/apps/HMMER/hmmalign";

 and the the current location of your vdjfasta installation:

  $self->{dnaVsegdb}    = "/tools/vdjfasta/db/imgt.VDJ.dna.nr.fa";
  $self->{dnaJsegdb}    = "/tools/vdjfasta/db/imgt.J.dna.fa";
  $self->{dnaDsegdb}    = "/tools/vdjfasta/db/imgt.HD.dna.nr.fa";
  $self->{dnaCsegdb}    = "/tools/vdjfasta/db/imgt.CH1.dna.nr.fa";

  $self->{VhVkHMM}      = "/tools/vdjfasta/db/Vh-linker-Vk.hmm";

Finally, add vdjfasta/bin to your path.

Description
===========

The scripts available in the package are all located in vdjfasta/bin. Terse usage descriptions are provided if the 
scripts are executed without arguments.

 bin: scripts available in the package
     fasta-vdj-pipeline.pl    complete analysis of antibody nucleotide sequences
     fasta-dna2ig.pl          identify and translate antibody coding frames
     fasta-getGermline.pl     identify segment composition
     fasta-scFv-align.pl      align antibody sequences
     fasta-vdj-sim.pl         perform antibody vdj simulations
     msa-restack.pl           perform CDR-restacking on scFv output alignment

Input
======

The fasta-vdj-pipeline.pl carries dna sequences from initial screening through full segment characterization, translation, and amino acid
identification of CDRs. An example of input and output is shown below:

Input would be a fasta file containing between 1 and 2,500 sequences. For larger datasets, a distributed solution is recommended.

>54_100
NNNNNNNNNNNNGGGCGATTGNTTTAGCGGCCGCGAATTCGCCCTTTGAAACACCTGTGGTTCTTCCTCCTCCTGGTGGC
AGCTCCCAGATGGGTCCTGTCTCAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGGGACCCTGTCCC
TCACCTGCGCTGTCTCTGGTGGCTCCATCAGCAGTAGTAACTGGTGGAGTTGGGTCCGCCAGCCCCCAGGGAAGGGGCTG
GAGTGGATTGGGGAAATCTATCATAGTGGGAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGA
CAAGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGAGAAAGC
TGGGGATTAAGTATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAGAGAGTCAGTCCTCCCCAACT
GTCTTCCCCCTCGTCTCCTGCGAGAGCCCCCTGTCTGATGAGAATTTGGTGGCCATGGGCTGCCTGGCCCGGGACTTCCT
GCCCAGCTCCATTTCCTTCTCCTGGAACTACCAGAACAACACTGAAGTCATGCAGGGTGTCAGAACCTTCCCAACACTGA
GGACAGGGGACAAATACACAGCTACCTCGCAGGTGTTACTGTCCGCCAAAAATGTCCTTGAAGGTTCAGATGAATACTTG
GTATGCAAAATCCACCATGGCAACAAAAATAAAGATCTGCATGTGCCGATTCCAGCTGTCGTTGAGATGAACCCCAATGT
GAGTGTGTTCATTCCACCACGTGATGCCTTCTCTGGCCCTGCACCCCGCAAGTCCAGACTCATCTGCGAGGCCACCAACT
TCAGTCCCAAACAGATCACAGTATCCTGGCTACAGGATGGGAAGCCTGTGAAATCTGGCTTCACNNCAGAGCCAGTGACT
GTCGANGCCNAANGNATCCAGACCCCAAANCTANNNGTCATNANCNNNCTGACCATCACTGAAAGNAGGGCNANTCNTTT
NANCNGCNGGACTAGTCCTTTANTGAGGGNNNNTGAGCTGNCGTANCATGNCATAGCTNNTCNGGNNNGNANTNNNNNTC
CNNNTCNNNANNNNNNANNNNNNNNNNNNNNNNNNNNNNNNCNGGGGNNCTANNNNNGNNCTNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNANNNTCNNNNNNNNNNNNN
>70_100
NNNNNNNNNNANNGGGCGATTGATTTAGCGGCCGCGAATTCGCCCTTTCCACGCTCCTGCTGCTGACCATCCCTTCATGG
GTCTTGTCCCAGATCACCTTGAAGGAGTCTGGTCCTACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCACCTT
CTCTGGGTTCTCACTCAGCACTAGTGGAGTGGGTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGGCCCTGGAGTGGCTTG
CACTCATTTATTGGAATGATGATAAGCGCTACAGCCCATCTCTGAAGAGCAGGCTCACCATCACCAAGGACACCTCCGAA
AACCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGACACAGCCACATATTACTGTGCACACGGATACAGCTATGG
TTACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAGAGAGTCAGTCCTCCCCAACTGTCTTCCCCCTCG
TCTCCTGCGAGAGCCCCCTGTCTGATGAGAATTTGGTGGCCATGGGCTGCCTGGCCCGGGACTTCCTGCCCAGCTCCATT
TCCTTCTCCTGGAACTACCAGAACAACACTGAAGTCATGCAGGGTGTCAGAACCTTCCCAACACTGAGGACAGGGGACAA
ATACACAGCTACCTCGCAGGTGTTACTGTCCGCCAAAAATGTCCTTGAAGGTTCAGATGAATACTTGGTATGCAAAATCC
ACCATGGCAACAAAAACAAAGATCTGCATGTGCCGATTCCAGCTGTCGTTGAGATGAACCCCAATGTGAGTGTGTTCATT
CCACCACGTGATGCCTTCTCTGGCCCTGCACCCCGCAAGTCCAGACTCATCTGCGAGGCCACCAACTTCAGTCCCAAACA
GATCACAGTATCCTGGCTACAGGATGGGAAGCCTGTGAAATCTGGCTTCACCACAGANCCAGTGACTGTCGANGCCAAAG
GATCCAGACCCCNAACCTACNNGGTCATAAGCACACTGACCATCACTGAAAGCANGGCNANNCGTTNNAANNTGCAGGAC
TAGTCCCTTTNNTGAGGTNATNNNNANCTNNNTANCATGNCATAGCTGTTNCTGGNNNGAAATGTNNCNGCTNNNNNNCN
NNNANNNNNNNNNNNNNNNNNNNNAAANNCNGGGNNGNCNANGNNNNNNTAANNCNNNNANNNGNNNNNNNNNNNNNNNN
NCNNNNNGGNNNNNNNNGCNNNNNNNCNNTAANGNANNNGGNNNNN

Output is a fasta file containing translated variable domains and headers describing the segment and CDR composition.

>54_100_frame_0;IGHV4-4 294 0;;IGHJ3 48 0;CARKLGIKYAFDIW;;IGHM 22 1
QVQLQESGPGLVKPSGTLSLTCAVSGGSISS.sNWWSWVRQPPGKGLEWIGEIYHSGSTNYNPSLKSRVTISVDKSKNQFSLKLSSVTAADTAVYYCARKLgi.kyaFDIWGQGTMVTVSS
>70_100_frame_2;IGHV2-5 297 1;IGHD5-5/5-18 18 0;IGHJ1 46 0;CAHGYSYGYFQHW;;IGHM 22 1
QITLKESGPTLVKPTQTLTLTCTFSGFSLSTsgVGVGWIRQPPGKALEWLALIYWNDDKRYSPSLKSRLTITKDTSENQVVLTMTNMDPVDTATYYCAHGYsy..gyFQHWGQGTLVTVSS

Validation
==========

Segment classification was validated against simulated datasets of somatically hypermutated segments of known composition. 


