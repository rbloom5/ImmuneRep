#!/usr/bin/python


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein


class Ab_read:
	def __init__(self, name = '', V = '', Vmut = [0], J='', Jmut=[0], ABtype = '',cdr3=Seq(''),cdr2=Seq(''),cdr1=Seq('')):
		self.name = name
		self.V = V
		self.Vmut = Vmut
		self.J = J
		self.Jmut = Jmut
		self.ABtype = ABtype
		self.cdr3 = cdr3
		self.cdr2 = cdr2
		self.cdr1 = cdr1

		if Vmut and Jmut:
			self.sh = Vmut[0]+Jmut[0]
		elif Vmut:
			self.sh = Vmut[0]
		else:
			self.sh = 0

	def __str__(self):
		return u'Read V={V}, J={J}, type = {type}'.format(
		V=self.V,
		J=self.J,
		type = self.ABtype
	)

	def __repr__(self):
		return self.__str__()



class Clone:
	def __init__(self, V = '', J='', cdr3=Seq(''), cdr2=Seq(''), cdr1=Seq(''), \
				ABtype = '', num_reads = None, percent_reads = None, IDs = [],\
				Vmut = None, Jmut = None, sh = None):
		self.V = V
		self.J = J
		self.cdr3 = cdr3
		self.cdr2 = cdr2
		self.cdr1 = cdr1
		self.num_reads = num_reads
		self.percent_reads = percent_reads
		self.ABtype = ABtype
		self.IDs = IDs
		self.Vmut = Vmut
		self.Jmut = Jmut
		self.sh = sh

	def __str__(self):
		return u'Clone V={V} J={J}, {num_reads} reads'.format(
			V=self.V,
			J=self.J,
			num_reads=self.num_reads,
		)

	def __repr__(self):
		return self.__str__()




class Cluster:
	def __init__(self, V = '', Js=[], cdr3s=[], cdr2s=[], cdr1s=[], \
				ABtypes = [], num_reads = None, percent_reads = None, \
				IDs = []):
		self.V = V
		self.Js = Js
		self.ABtypes = ABtypes
		self.cdr3s = cdr3s
		self.cdr2s = cdr2s
		self.cdr1s = cdr1s
		self.num_reads = num_reads
		self.percent_reads = percent_reads
		self.IDs = IDs

	def __str__(self):
		return u'Cluster V={V}, {num_reads} total reads'.format(
		V=self.V,
		num_reads=self.num_reads,
	)

	def __repr__(self):
		return self.__str__()