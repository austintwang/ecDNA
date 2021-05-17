def separate_repeats(seqs):
	reads_repeat = {}
	reads_nonrepeat = {}
	reads_foldback = {}

	for k, v in seqs.items():
		has_repeat = False
		has_foldback = False
		has_nonrepeat = False
		for i in range(1, len(v)):
			if v[i][3] != v[i-1][3]:
				has_nonrepeat = True
				continue
			overlap = 
			if 

