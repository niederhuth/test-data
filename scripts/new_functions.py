import gzip
import gffutils
import pandas as pd
import pybedtools as pbt
from os import remove
from numpy import float64
from itertools import product

#Read allc file and convert to bedfile for use with bedtools
def allc2bed(allc, bedfile):
	#Check if file is compressed with gzip and open
	if allc.endswith('gz'):
		#Open gzip file as text
		in_file = gzip.open(allc, 'rt')
	else:
		in_file = open(allc, 'r')
	#Open output bedfile
	out_file = open(bedfile, 'w')
	#Read in each line
	for line in in_file.readlines():
		#Split the line on tab
		a = line.split('\t')
		#Add additional fields to make it a bed-file
		out_file.write('\t'.join([a[0], a[1], a[1], '.', '.'] + a[2:7]))
	#Close open file handles
	in_file.close()
	out_file.close()

#Interpret sequence context, taken from methylpy utilities
def expand_nucleotide_code(mc_type=['C']):
	iub_dict = {'N':['A','C','G','T','N'],
				'H':['A','C','T','H'],
				'D':['A','G','T','D'],
				'B':['C','G','T','B'],
				'A':['A','C','G','A'],
				'R':['A','G','R'],
				'Y':['C','T','Y'],
				'K':['G','T','K'],
				'M':['A','C','M'],
				'S':['G','C','S'],
				'W':['A','T','W'],
				'C':['C'],'G':['G'],'T':['T'],'A':['A']}
	mc_class = list(mc_type) # copy
	if 'C' in mc_type:
		mc_class.extend(['CGN', 'CHG', 'CHH','CNN'])
	elif 'CG' in mc_type:
		mc_class.extend(['CGN'])
	mc_class_final = []
	for motif in mc_class:
		mc_class_final.extend([''.join(i) for i in product(*[iub_dict[nuc] for nuc in motif])])
	return(set(mc_class_final))

#function for filtering annotation files based on feature (gene, exon, mRNA, etc)
def feature_filter(x,feature):
	if feature:
		return x[2] == feature
	else:
		return x

#function for filtering annotation files based on chromosome
def chr_filter(x,chr):
	return x.chrom in chr

#Collect mC data for a context
#If 'site_cutoff_only=True', then the cutoff will only apply to 
#tallying the number of sites & methylated reads called by methylpy
def get_mC_data(df, mc_type='C', cutoff=0, site_cutoff_only=False):
	mC1 = mC2 = mC3 = mC4 = 0
	#expand nucleotide list for a given context
	exp_mc_type = expand_nucleotide_code(mc_type=[mc_type])
	#Filter input pandas dataframe for correct mc_class
	filt_df = df[df['mc_class'].isin(exp_mc_type)]
	#Count number of sites passing the cutoff
	mC1 = len(filt_df[filt_df['total'] > cutoff])
	#if no sites for that context, set all values to 'NA'
	if mC1 == 0:
		mC1 = mC2 = mC3 = mC4 = 'NA'
	else:
		#Count number of sites called as methylated (column 7 of allc) by methylpy
		mC2 = sum(filt_df[filt_df['total'] > cutoff]['methylated'])
		#Check if site_cutoff_only is set to true or false
		if site_cutoff_only:
			mC3 = sum(filt_df['mc_count'])
			mC4 = sum(filt_df['mc_count'])
		else:
			mC3 = sum(filt_df[filt_df['total'] > cutoff]['mc_count'])
			mC4 = sum(filt_df[filt_df['total'] > cutoff]['mc_count'])
		#if no sites for that context, set to 'NA'
		if mC1 == 0:
			mC1 = mC2 = mC3 = mC4 = 'NA'
	#Return results
	return([mc_type,mC1,mC2,mC3,mC4])

#Calculate methylation levels for features. These can be genes, transposons, etc.
#'mC_bed' is a bedfile created from allc2bed
#'gff' is a gff or gtf file
#'mc_type' are the specific sequence contexts to look at. Supports expanded nucleotide code, 
#e.g. 'CHH'
#'primary_feature' is the main feature to examine, e.g. 'gene'
#'secondary_feature' is a second sub-feature used to filter data to specific parts of the primary 
#feature. For example, if you want to restrict analysis of genic methylation to coding sequences 
#(CDS) or exons.
#
def feature_methylation(mC_bed, gff, output, mc_type=['CG','CHG','CHH'], primary_feature='gene', secondary_feature=(), cutoff=0, site_cutoff_only=False):
	out_file = open(output, 'w')
	columns = ['Total_C', 'Methylated_C', 'Total_Reads', 'Methylated_Reads', 'Weighted_mC']
	header = ['Feature']
	for i in mc_type:
		for i2 in columns:
			header = header + [i + '_' + i2]
	out_file.write('\t'.join(header) + '\n')
	db = gffutils.create_db(gff, dbfn='tmp.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
	tmp = pbt.BedTool(mC_bed)
	for pf in db.features_of_type(primary_feature):
		if secondary_feature:
			sub_feat = ""
			for i in db.children(pf, featuretype=secondary_feature, order_by='start'):
				sub_feat = sub_feat + str(i) + '\n'
			pf_bed = pbt.BedTool(sub_feat, from_string=True)
		else:
			pf_bed = pbt.BedTool(pf, from_string=True)
		mapping = pd.read_csv(pbt.bedtool.BedTool.intersect(tmp, pf_bed).fn, usecols=[6,7,8,9], sep="\t", names=['mc_class','mc_count','total','methylated'], dtype={'mc_count':int,'total':int,'methylated':int})	
		out_row = [pf.id]
		for i in mc_type:
			mC = get_mC_data(mapping, mc_type=i, cutoff=cutoff, site_cutoff_only=site_cutoff_only)
			if mC[1] == 'NA':
				out_row += ['NA','NA','NA','NA','NA']
			else:
				out_row += [mC[1], mC[2], mC[3], mC[4], (float64(mC[4])/float64(mC[3]))]
		out_file.write('\t'.join([str(elem) for elem in out_row]) + '\n')
	remove('tmp.db')

def feature_methylation(mC_bed, gff, output, mc_type=['CG','CHG','CHH'], primary_feature='gene', secondary_feature=(), cutoff=0, site_cutoff_only=False):
	out_file = open(output, 'w')
	columns = ['Total_C', 'Methylated_C', 'Total_Reads', 'Methylated_Reads', 'Weighted_mC']
	header = ['Feature']
	for i in mc_type:
		for i2 in columns:
			header = header + [i + '_' + i2]
	out_file.write('\t'.join(header) + '\n')
	db = gffutils.create_db(gff, dbfn='temp.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
	tmp = pbt.BedTool(mC_bed)
	tmp2 = pbt.BedTool(gff).filter(feature_filter,secondary_feature).filter(chr_filter,chrs)
	mapping = pbt.bedtool.BedTool.intersect(tmp, tmp2, wa=True, wb=True)
	for pf in db.features_of_type(primary_feature):
		

	feat = ""
	for pf in db.features_of_type(primary_feature):
		if secondary_feature:
			for i in db.children(pf, featuretype=secondary_feature, order_by='start'):
				feat = feat + str(i) + '\n'
		else:
			feat = feat + str(pf) + '\n'
	pf_bed = pbt.BedTool(feat, from_string=True)
	tmp = pbt.BedTool(mC_bed)
	mapping = pbt.bedtool.BedTool.intersect(tmp, pf_bed, wa=True, wb=True)
	a = pd.read_csv(mapping.fn, usecols=[6,7,8,9,20], sep="\t", names=['mc_class','mc_count','total','methylated'], dtype={'mc_count':int,'total':int,'methylated':int})	
	


		out_row = [pf.id]
		for i in mc_type:
			mC = get_mC_data(mapping, mc_type=i, cutoff=cutoff, site_cutoff_only=site_cutoff_only)
			if mC[1] == 'NA':
				out_row += ['NA','NA','NA','NA','NA']
			else:
				out_row += [mC[1], mC[2], mC[3], mC[4], (float64(mC[4])/float64(mC[3]))]
		out_file.write('\t'.join([str(elem) for elem in out_row]) + '\n')
	remove('tmp.db')
