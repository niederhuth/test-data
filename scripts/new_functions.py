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
#This function is essential for interpreting methylation contexts
#that include non-standard letters, such as 'CHH'
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
	mc_class = list(mc_type)
	if 'C' in mc_type:
		mc_class.extend(['CNN'])
	elif 'CG' in mc_type:
		mc_class.extend(['CGN'])
	elif 'CH' in mc_type:
		mc_class.extend(['CHN'])
	mc_class_final = []
	for motif in mc_class:
		mc_class_final.extend([''.join(i) for i in product(*[iub_dict[nuc] for nuc in motif])])
	return(set(mc_class_final))

#Filter annotation files based on feature (gene, exon, etc)
#Largely intended for use with pybedtools
#'x' is a line from a gff file
#'feature' is the particular feature type to keep (gene, exon, etc)
def feature_filter(x,feature):
	if feature:
		return x[2] == feature
	else:
		return x

#Filter files, keeping those in the chromosomes listed
#Largely intended for use with pybedtools
#'x' is a line from the file, e.g. a gff or bed file
#'chr' can be a list of chromosome names. These must match those
#in the record
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
			mC3 = sum(filt_df['total'])
			mC4 = sum(filt_df['mc_count'])
		else:
			mC3 = sum(filt_df[filt_df['total'] > cutoff]['total'])
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
#'cutoff' applies a minimum read cutoff to the data. Any site below that cutoff is not counted.
#'site_cutoff_only', when set to True, will only apply the cutoff to the summing up the number of
#sites, but will calculated the total weighted methylation using all the data regardless of the cutoff.
#'chrs' is a list of chromosomes or sequences (e.g. scaffolds/contigs) to keep.
def feature_methylation(mC_bed, gff, output, mc_type=['CG','CHG','CHH'], primary_feature='gene', 
	secondary_feature=(), chrs=[], cutoff=0, site_cutoff_only=False):
	#Open output file for reading
	out_file = open(output, 'w')
	#Build the header row
	header = ['Feature']
	#List of output columns to be appended to each mc_type
	columns = ['Total_C', 'Methylated_C', 'Total_Reads', 'Methylated_Reads', 'Weighted_mC']
	#Iterate over each mc_type
	for i in mc_type:
		#Iterate over each of the columns
		for i2 in columns:
			#Combine the mc_type with the column and add this to the header row
			header = header + [i + '_' + i2]
	#Write the header row to the output file
	out_file.write('\t'.join(header) + '\n')
	#Create a temporary database of the gff file using gffutils. This will allow us to easily
	#map child features, like CDS or exon, to a parent feature. 
	db = gffutils.create_db(gff, dbfn='temp.db', force=True, keep_order=True, 
		merge_strategy='merge', sort_attribute_values=True)
	#Create an empty variable called sf_bed (secondary feature bed), from which we will build
	#a new bed file, of secondary feature coordinates, but with the primary feature id as name
	sf_bed = ''
	#Create an empty list that will be populated with the primary feature ids
	pf_list = []
	#Iterate over each primary feature in the gff database
	for pf in db.features_of_type(primary_feature):
		#Check if the primary feature is located on the chromosomes being analyzed,
		#otherwise that primary feature is skipped.
		if pf.seqid in chrs:
			#Add the primary feature id to our list
			pf_list = pf_list + [pf.id]
			#If a secondary feature is provided, we will use the coordinates from these to map the 
			#methylation data. This is so we can exclude unwanted regions of the primary feature, 
			#such as intronic sequences. However, we still want to use the primary feature name,
			#so we create a new bed file, with this name instead.
			if secondary_feature:
				#Loop over each secondary feature that is a child of the primary feature
				for i in db.children(pf, featuretype=secondary_feature, order_by='start'):
					#Create our bed file, by joining together the sequence id (chromosome), start, stop,
					#and primary feature id
					sf_bed = sf_bed + '\t'.join([str(i.seqid), str(i.start), str(i.stop), str(pf.id)]) + '\n'
			#if a secondary feature is not provided, We again create a new bed file as above, but using,
			#using the coordinates of the primary feature. 
			else:
				sf_bed = sf_bed + '\t'.join([str(pf.seqid), str(pf.start), str(pf.stop), str(pf.id)]) + '\n'
	#We then create pybedtools objects for both our input mC_bed from allc2bed and our sf_bed.
	#We then find the intersection of these two bed file, outputing the records for both
	map_bed = pbt.bedtool.BedTool.intersect(pbt.BedTool(mC_bed), 
		pbt.BedTool(sf_bed, from_string=True), wa=True, wb=True)
	#We now read in the mapped methylation data (map_bed) as a pandas frame, keeping only the relevant 
	#information to reduce the size of our dataframe. Specifically we keep column 6 (mc_class), 
	#column 7 (count of methylated reads), column 8 (total read count), column 9 (whether site is called
	#methylated '1' or unmethylated '0' by methylpy), and column 13, our primary feature id
	map_df = pd.read_csv(map_bed.fn, usecols=[6,7,8,9,13], sep="\t", 
		names=['mc_class','mc_count','total','methylated','id'], 
		dtype={'mc_count':int,'total':int,'methylated':int,'id':str})
	#We now loop over each primary feature id in our pf_list
	for pf in pf_list:
		#We extract the rows of mapped methylation data for that primary feature id from
		#from the pandas dataframe
		feat_map = map_df[map_df['id']== pf]
		#We initiate a new list called 'out_row', the first element being the primary feature id
		out_row = [pf]
		#For each mc_type, we then add up the methylation data for the primary feature
		for i in mc_type:
			#We use the function 'get_mC_data' for this task
			mC = get_mC_data(feat_map, mc_type=i, cutoff=cutoff, site_cutoff_only=site_cutoff_only)
			#If for that mc_type, there is not data in the total reads column, we set all the other
			#values to "NA" and add it to our 'out_row' list.
			if mC[1] == 'NA':
				out_row += ['NA','NA','NA','NA','NA']
			#Otherwise, we add the methylation data to the 'out_row' list, and also calculate
			#the weighted methylation by dividing the methylated reads column by the total reads
			#column. We using numpy.float64, to avoid certain errors that can arise
			else:
				out_row += [mC[1], mC[2], mC[3], mC[4], (float64(mC[4])/float64(mC[3]))]
		#We then write the data to our output file
		out_file.write('\t'.join([str(elem) for elem in out_row]) + '\n')
	#When we are done, we close the output file
	out_file.close()
	#And we delete our temporary gff database.
	remove('temp.db')

