#!/usr/bin/python




class Rep_stats:

	def __init__(self, file_locs):
		self.features={}

		for f in file_locs:
			f_id = f.split('/')[-1][:-4]
			with open(f) as rep:
				cur_rep = pickle.load(rep)
				self.features[f_id] = cur_rep.features

	

