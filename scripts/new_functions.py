import os
import sys
import gzip
import itertools

import numpy as np
import pandas as pd
import pybedtools as pbt

#interpret sequence context, taken from methylpy utilities
def expand_nucleotide_code(mc_type=['C']):
	iub_dict = {'N':['A','C','G','T','N'],'H':['A','C','T','H'],'D':['A','G','T','D'],'B':['C','G','T','B'],'A':['A','C','G','A'],'R':['A','G','R'],'Y':['C','T','Y'],'K':['G','T','K'],'M':['A','C','M'],'S':['G','C','S'],'W':['A','T','W'],'C':['C'],'G':['G'],'T':['T'],'A':['A']}
	mc_class = list(mc_type) # copy
	if 'C' in mc_type:
		mc_class.extend(['CGN', 'CHG', 'CHH','CNN','CH','CN'])
	elif 'CG' in mc_type:
		mc_class.extend(['CGN'])
	mc_class_final = []
	for motif in mc_class:
		mc_class_final.extend([''.join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in motif])])
	return(set(mc_class_final))

#function for filtering annotation files based on feature (gene, exon, mRNA, etc)
def feature_filter(x,feature):
	if feature:
		return x[2] == feature
	else:
		return x

#function for filtering annotation files based on strand
def strand_filter(x,strand):
	return x.strand == strand

#function for filtering annotation files based on chromosome
def chr_filter(x,chr):
	return x.chrom not in chr

#Read allc file and convert to bedfile
def allc2bed(allc):
	#get first line
	if allc.endswith('gz'):
		header = gzip.open(allc).readline().rstrip()
	else:
		header = open(allc).readline().rstrip()
	#check if first line contains header
	if header == 'chr\tpos\tstrand\tmc_class\tmc_count\ttotal\tmethylated':
		#read in allc file to pandas dataframe
		a = pd.read_table(allc,dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int})
	else:
		#read in allc file to pandas dataframe, add header if does not have one
		a = pd.read_table(allc,names=['chr','pos','strand','mc_class','mc_count','total','methylated'],dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int})
	#add new columns
	a['pos2'] = a['pos']
	a['name'] = a.index
	a['score'] = '.'
	#reorder columns
	a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
	#return bed file
	a = pbt.BedTool.from_dataframe(a)
	return a

#Collect mC data for a context
def get_mC_data(a,mc_type='C',cutoff=0):
	#expand nucleotide list for a given context
	b = expand_nucleotide_code(mc_type=[mc_type])
	d1 = d2 = d3 = d4 = 0
	#iterate over each line
	for c in a.itertuples():
		#check if line is correct context
		if c[4] in b:
			#check if meets minimum cutoff for read depth
			if int(c[6]) >= int(cutoff):
				#add up number of sites
				d1 = d1 + 1
				#add up number of sites called methylated by methylpy
				d2 = d2 + int(c[7])
				#add up total reads covering a site
				d3 = d3 + int(c[6])
				#add up total methylated reads covering a site
				d4 = d4 + int(c[5])
	#create list
	e = [mc_type,d1,d2,d3,d4]
	#return that list
	return e

#for calculating methylation levels in windows across the genome
def genome_window_methylation(allc,genome_file,output=(),mc_type=['CG','CHG','CHH'],window_size=100000,stepsize=50000,cutoff=0,chrs=[]):
	#read in allc file
	a = allc2bed(allc)
	#create output data frame
	c = ['Chr','Window']
	columns=['Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#make windows
	w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file).filter(chr_filter,chrs),g=genome_file,w=window_size,s=stepsize,i='srcwinnum')
	#intersect bedfiles with pybedtools
	mapping = pbt.bedtool.BedTool.intersect(a,w_bed,wa=True,wb=True)
	del(w_bed,a)
	#convert to pandas dataframe
	m = pd.read_table(mapping.fn,header=None,usecols=[13,6,7,8,9])
	del(mapping)
	#split srcwinnum
	b = m[13].str.split('_', n = 1, expand = True)
	#make new columns from srcwinnum
	m['Chr'] = f[0]
	m['Window'] = f[1]
	del(f)
	#reorder data frame
	m = m[['Chr','Window',13,6,7,8,9]]
	#iterate over each chromosome
	for g in chrs:
		#get windows for that specific chromosome
		windows = list(m[m['Chr'].isin([str(g)])]['Window'].drop_duplicates())
		#iterate over each window in each chromosome
		for h in windows:
			#filter for rows matching specific chr & window number
			i = m[m['Chr'].isin([str(g)]) & m['Window'].isin([str(h)])]
			#make list for methylation data
			j = [g,h]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				l = get_mC_data(i,mc_type=k,cutoff=cutoff)
				#delete first line of list
				del(l[0])
				#Calculate weighted methylation and add this to list of data for other mc_types
				j = j + l + [(np.float64(l[4])/np.float64(l[3]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b

def metaplot(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],window_number=60,updown_stream=2000,feature=(),cutoff=0,chrs=[]):
	#read in allc file
	a = allc2bed(allc)
	#create output data frame
	c = ['Window']
	columns=['Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#read annotation file and filter on feature
	f_bed = pbt.BedTool(annotations).filter(feature_filter,feature).filter(chr_filter,chrs).saveas('f_bed.tmp')
	#if updown_stream set to 0, only region to look at is the feature itself, e.g. f_bed
	if updown_stream == 0:
		regions=[f_bed]
	#if updown_stream specified, create bed files for upstream regions (u_bed) and down stream regions (d_bed)
	else:
		u_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=updown_stream,r=0,s=True).saveas('u_bed.tmp')
		d_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=0,r=updown_stream,s=True).saveas('d_bed.tmp')
		regions=[u_bed,f_bed,d_bed]
	#set window number to 1
	window = 1
	#iterate over each region and collect methylation data
	for f in regions:
		#filter bed files based on strand
		p_bed = f.filter(strand_filter,strand='+').saveas('p_bed.tmp')
		n_bed = f.filter(strand_filter,strand='-').saveas('n_bed.tmp')
		#make windows
		pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(window_number/3),i='srcwinnum')
		nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,n=int(window_number/3),i='srcwinnum',reverse=True)
		#combine windows
		w_bed = pw_bed.cat(nw_bed, postmerge=False)
		#intersect bedfiles with pybedtools
		mapping = pbt.bedtool.BedTool.intersect(a,w_bed,wa=True,wb=True)
		del(w_bed,pw_bed,nw_bed)
		#convert to pandas dataframe
		m = pd.read_table(mapping.fn,header=None,usecols=[13,6,7,8,9])
		del(mapping)
		#split srcwinnum
		g = m[13].str.split('_', n = 1, expand = True)
		#make new columns from srcwinnum
		m['Name'] = g[0]
		m['Window'] = g[1]
		del(g)
		#reorder data frame
		m = m[['Name','Window',13,6,7,8,9]]
		#iterate list of window numbers
		for h in list(range(1,int(window_number/3)+1)):
			#filter for rows matching specific window number
			i = m[m['Window'].isin([str(h)])]
			#make list for methylation data
			j = [window]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				l = get_mC_data(i,mc_type=k,cutoff=cutoff)
				#Calculate weighted methylation and add this to list of data for other mc_types
				j = j + [(np.float64(l[4])/np.float64(l[3]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
			#count windows
			window = window + 1
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
	#remove temporary files created
	tmp=['f_bed.tmp','u_bed.tmp','d_bed.tmp','p_bed.tmp','n_bed.tmp']
	for n in tmp:
		os.remove(n)
