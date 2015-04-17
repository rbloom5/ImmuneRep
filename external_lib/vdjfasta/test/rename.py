from Bio import SeqIO
#from StringIO import StringIO

#handle = StringIO("testing.fasta")
def rename():
	SeqIO.convert("SRR1298742.fastq","fastq",handle,"fasta")
	print "5"

rename()