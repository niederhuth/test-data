import gzip
import gffutils
import pandas as pd
import pybedtools as pbt
from itertools import product

#Read allc file and convert to bedfile for use with bedtools
def allc2bed(allc,bedfile):
	#Check if file is compressed with gzip and open
	if allc.endswith('gz'):
		#Open gzip file as text
		in_file = gzip.open(allc, 'rt')
	else:
		in_file = open(allc, 'r')
	#Open output bedfile
	out_file=open(bedfile,'w')
	#Read in each line
	for line in in_file.readlines():
		#Split the line on tab
		a = line.split('\t')
		#Output the file as tab delimited with additional elements to
		#make it a bedfile 
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

#Collect mC data for a context
#If 'site_cutoff_only=True', then the cutoff will only apply to 
#tallying the number of sites & methylated reads called by methylpy
def get_mC_data(df,mc_type='C',cutoff=0,site_cutoff_only=False):
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
#'mc_type' are the specific sequence contexts to look at. Supports expanded nucleotide code, e.g. 'CHH'
#'primary_feature' is the main feature to examine, e.g. 'gene'
#'secondary_feature' is a second sub-feature used to filter data to specific parts of the primary feature.
#for example, if you want to restrict analysis of genic methylation to coding sequences (CDS) or exons.
#
def feature_methylation(mC_bed,gff,mc_type=['CG','CHG','CHH'],primary_feature='gene',secondary_feature=()):
	#Create a gff database with gffutils
	db = gffutils.create_db(gff, dbfn='temp.db', force=True, keep_order=True, 
		merge_strategy='merge', sort_attribute_values=True)
	#Iterate over each primary feature in the database
	for pf in db.features_of_type(primary_feature):
		if secondary_feature:
			for i in db.children(gene, featuretype='cds', order_by='start'):
		#Check if a secondary_feature is defined, if not, default to only using the primary feature
		if secondary_feature:
			for i in db.children(gene, featuretype='cds', order_by='start'):

		else:
			#Create a pybedtools object from the primary feature
			pf_bed = pbt.BedTool(pf,from_string=True)
		for i in db.children(gene, featuretype='cds', order_by='start'):
			for i in db.children(gene, featuretype='cds', order_by='start'):
		b = pbt.BedTool(a,from_string=True)
		m = pd.read_csv(pbt.bedtool.BedTool.intersect(allc,b).fn,usecols=[6,7,8,9],sep="\t",
			names=['mc_class','mc_count','total','methylated'],
			dtype={'mc_count':int,'total':int,'methylated':int})
		for x in mc_type:
			print(get_mC_data(m,mc_type=x,cutoff=0,site_cutoff_only=False)
		a=""
