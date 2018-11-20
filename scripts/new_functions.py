#Collect mC data for a context
def get_mC_data(a,mc_type=['C'],cutoff=0):
	b = expand_nucleotide_code(mc_type)
	d1 = d2 = d3 = d4 = 0
	for c in a.itertuples():
		if c[4] in b:
			if int(c[6]) >= int(cutoff):
				d1 = d1 + 1
				d2 = d2 + int(c[7])
				d3 = d3 + int(c[6])
				d4 = d4 + int(c[5])
	e = [mc_type,d1,d2,d3,d4]
	return e

#
def something():
