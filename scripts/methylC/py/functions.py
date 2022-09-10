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


#Filter annotation files based on strand
#Largely intended for use with pybedtools
#'x' is a line from a bed file
#'strand' is the strand you want to keep, either '+' or '-'
def strand_filter(x,strand):
	return x.strand == strand

#Filter files based on feature names, e.g. gene names
#Largely intended for use with pybedtools
#'x' is a line from a bed file
#name_list is a python list of feature names to keep
def name_filter(x,name_list):
	return x.name in name_list


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


#This function creates a bed file where one feature (called the primary_feature or pf) has been
#filtered to only include regions from another sub-feature (called the secondary_feature or sf).
#For example, if provided a gff file and the the primary_feature is "gene" and the secondary_feature 
#is "exon", then this script will create a set of bed coordinates where each exon (secondary_feature)
#has the the name of the corresponding gene (primary_feature). These are used mostly in the map_to_gff
#and metaplot functions, when we want methylation levels that only correspond to certain parts of a 
#gene (e.g. exon, CDS, intron). We list it as a separate function so we can call it in other functions.
#'gff' is a gff or gtf file.
#'primary_feature' is the main feature to examine, e.g. 'gene'
#'secondary_feature' is a second sub-feature used to filter data to specific parts of the primary 
#feature. For example, if you want to restrict analysis of genic methylation to coding sequences 
#(CDS) or exons.
#'chrs' is a list of chromosomes or sequences (e.g. scaffolds/contigs) to keep.
#return_pf_bed, if set to TRUE, will also return a bed file of the primary_feature
#This function always returns the result and does not output to file.
def sf2pf_bed(gff, primary_feature='gene', secondary_feature=(), chrs=[], return_pf_bed=False):
	#Create a temporary database of the gff file using gffutils. This will allow us to easily
	#map child features, like CDS or exon, to a parent feature. 
	db = gffutils.create_db(gff, dbfn='temp.db', force=True, keep_order=True, 
		merge_strategy='merge', sort_attribute_values=True)
	#Create an empty variable called sf_bed (secondary feature bed), from which we will build
	#a new bed file, of secondary feature coordinates, but with the primary feature id as name
	sf_bed = pf_bed = ''
	#Create an empty list that will be populated with the primary feature ids
	pf_list = []
	#Iterate over each primary feature in the gff database
	for pf in db.features_of_type(primary_feature):
		#Check if the primary feature is located on the chromosomes being analyzed,
		#otherwise that primary feature is skipped.
		if pf.seqid in chrs:
			#Add the primary feature id to our list
			pf_list = pf_list + [pf.id]
			#Create a bed file of primary of the primary feature
			pf_bed = pf_bed + '\t'.join([str(pf.seqid), str(pf.start-1), str(pf.stop), str(pf.id), 
				str(pf.score), str(pf.strand)]) + '\n'
			#If a secondary feature is provided, we will use the coordinates from these to map the 
			#methylation data. This is so we can exclude unwanted regions of the primary feature, 
			#such as intronic sequences. However, we still want to use the primary feature name,
			#so we create a new bed file, with this name instead.
			if secondary_feature:
				#Loop over each secondary feature that is a child of the primary feature
				for sf in db.children(pf, featuretype=secondary_feature, order_by='start'):
					#Create our bed file, by joining together the sequence id (chromosome), start, stop,
					#and primary feature id
					sf_bed = sf_bed + '\t'.join([str(sf.seqid), str(sf.start-1), str(sf.stop), str(pf.id),
						str(sf.score), str(sf.strand)]) + '\n'
	#Read the sf_bed/pf_bed with pybedtools, saving a temporary version of this, to text.
	#We save a temporary version of this, because sometimes, the bedfile iterator from pybedtools will
	#get deleted.
	#if secondary_feature is provided, we read in sf_bed
	if secondary_feature:
		sf_pbt_bed = pbt.BedTool(sf_bed, from_string=True).saveas('sf_bed.tmp')
	#if a secondary feature is not provided, We again create a new bed file as above, but using,
	#using the coordinates of the pf_bed.
	else:
		sf_pbt_bed = pbt.BedTool(pf_bed, from_string=True).saveas('sf_bed.tmp')
	#Check if return_pf_bed is TRUE
	if return_pf_bed:
		#if flank_distance or return_bf_bed, then create a pybedtools object from that.
		pf_pbt_bed = pbt.BedTool(pf_bed, from_string=True).saveas('pf_bed.tmp')
		return(sf_pbt_bed,pf_list,pf_pbt_bed)
	else:
		return(sf_pbt_bed,pf_list)


#Calculate methylation levels for features in a bed file. 
#These can be anything, as long as the feature name is in the 4th column of the bed file.
#'mC_bed' is a bedfile created from allc2bed
#'feat_bed' is a bedfile of the features you want to map to
#'feature_ids' is an optional list of specific feature ids you want to look at. It needs to be a 
#python list, e.g. in the format ['x','y','z']. It can be used to restict analyses to 
#a subset of features, or included in case data may be missing, which will guarantee 
#that these features are always in the output file, even with an 'NA' value. 
#'mc_type' are the specific sequence contexts to look at. Supports expanded nucleotide code, 
#e.g. 'CHH'
#'cutoff' applies a minimum read cutoff to the data. Any site below that cutoff is not counted.
#'site_cutoff_only', when set to True, will only apply the cutoff to the summing up the number of
#sites, but will calculated the total weighted methylation using all the data regardless of the cutoff.
#'chrs' is a list of chromosomes or sequences (e.g. scaffolds/contigs) to keep.
def map_to_bed(mC_bed, feat_bed, output, feature_ids=[], mc_type=['CG','CHG','CHH'], 
	cutoff=0, site_cutoff_only=False):
	#If a feature_ids is not provided, we build one using column 4 of the feat_bed.
	if not feature_ids:
		feature_ids = list(pd.unique(pd.read_csv(pbt.BedTool(feat_bed).fn, sep='\t', 
			header=None)[3]).astype(str))
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
	#We create pybedtools objects for both our input mC_bed  and our feat_bed.
	#We then find the intersection of these two bed filed, outputing the records for both
	map_bed = pbt.bedtool.BedTool.intersect(pbt.BedTool(mC_bed), 
		pbt.BedTool(feat_bed), wa=True, wb=True)
	#We now read in the mapped methylation data (map_bed) as a pandas frame, keeping only the relevant 
	#information to reduce the size of our dataframe. Specifically we keep column 6 (mc_class),
	#column 7 (count of methylated reads), column 8 (total read count), column 9 (whether site is called
	#methylated '1' or unmethylated '0' by methylpy), and column 13, our primary feature id
	map_df = pd.read_csv(map_bed.fn, usecols=[6,7,8,9,13], sep="\t", 
		names=['mc_class','mc_count','total','methylated','id'], 
		dtype={'mc_count':int,'total':int,'methylated':int,'id':str})
	#We now loop over each primary feature id in our pf_list
	for feat in feature_ids:
		#We extract the rows of mapped methylation data for that primary feature id from
		#from the pandas dataframe
		feat_map = map_df[map_df['id'] == feat]
		#We initiate a new list called 'out_row', the first element being the primary feature id
		out_row = [feat]
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


#Calculate methylation levels for features in a gff file. 
#This command is an extension of the 'map_to_bed' function, for gff files.
#It also allows one to filter data first by a secondary sub-feature, e.g. 'exon', before
#mapping to a primary feature, e.g. 'gene'
#'mC_bed' is a bedfile created from allc2bed
#'gff' is a gff or gtf file
#'output' is the output file
#'mc_type' are the specific sequence contexts to look at. Supports expanded nucleotide code, 
#e.g. 'CHH'
#'primary_feature' is the main feature to examine, e.g. 'gene'
#'secondary_feature' is a second sub-feature used to filter data to specific parts of the primary 
#feature. For example, if you want to restrict analysis of genic methylation to coding sequences 
#(CDS) or exons.
#'chrs' is a list of chromosomes or sequences (e.g. scaffolds/contigs) to keep.
#'feature_ids' is an optional list of specific feature ids you want to look at. It needs to be a 
#python list, e.g. in the format ['x','y','z']. It can be used to restict analyses to 
#a subset of features, or included in case data may be missing, which will guarantee 
#that these features are always in the output file, even with an 'NA' value. 
#'cutoff' applies a minimum read cutoff to the data. Any site below that cutoff is not counted.
#'site_cutoff_only', when set to True, will only apply the cutoff to the summing up the number of
#sites, but will calculated the total weighted methylation using all the data regardless of the cutoff.
def map_to_gff(mC_bed, gff, output, mc_type=['CG','CHG','CHH'], primary_feature='gene', 
	secondary_feature=(), chrs=[], feature_ids=[], cutoff=0, site_cutoff_only=False):
	#Use sf2pf_bed to map secondary_feature to primary_feature and create
	sf_pbt_bed,pf_list = sf2pf_bed(gff, primary_feature=primary_feature, 
		secondary_feature=secondary_feature, chrs=chrs, return_pf_bed=False)
	#If feature_ids are provided, then all analyses will be restricted to that list of ids
	if feature_ids:
		#Use 'name_filter' function to filter sf_pbt_bed for those features
		sf_pbt_bed = sf_pbt_bed.filter(name_filter,feature_ids)
	#We now invoke map_to_bed to map mC data to our sf_pbt_bed and call methylation for each feature
	#in our pf_list. This will output to our 'output' file.
	map_to_bed(mC_bed, sf_pbt_bed, output, feature_ids=pf_list, mc_type=mc_type, 
		cutoff=cutoff, site_cutoff_only=site_cutoff_only)
	#And we delete our temporary gff database and bed file.
	for tmp_file in ['temp.db','temp.bed']:
		remove(tmp_file)


#Calculate methylation levels going across a set of genomic features (e.g. genes, transposons) using a set
#number of windows. This will average the weighted methylation across all those features in those windows. 
#So if you look at all genes, each window will have data from all genes for that specific region. 
#Note that this function does not plot data as a set number of basepairs (e.g. from 0-500 bp, etc). 
#In other words, a window for the gene body may be variable from gene to gene depending on the size of 
#the gene.
#'mC_bed' is a bedfile created from allc2bed
#'gff' is a gff or gtf file
#genome_file is either a fai (generated by samtools faidx) or .genome file used to get a list of chromosomes)
#'output' is the output file
#'mc_type' are the specific sequence contexts to look at. Supports expanded nucleotide code, 
#e.g. 'CHH'
#'flank_distance' is the number of basepairs upstream or downstream to include in the plot. If set to 
#zero, then only the primary_feature will be analyzed
#'window_number' is how many windows you want to use. This number will be applied to each region analysed. 
#In other words, if 'window_number' = 20, then the upstream region will be made up of 20 windows, the 
#primary feature will be broken down to 20 windows, the downstream region will be broken down to 20 windows
#'primary_feature' is the main feature to examine, e.g. 'gene'
#'secondary_feature' is a second sub-feature used to filter data to specific parts of the primary 
#feature. For example, if you want to restrict analysis of genic methylation to coding sequences 
#(CDS) or exons.
#'chrs' is a list of chromosomes or sequences (e.g. scaffolds/contigs) to keep.
#'feature_ids' is an optional list of specific feature ids you want to look at. It needs to be a 
#python list, e.g. in the format ['x','y','z']. It can be used to restict analyses to 
#a subset of features, or included in case data may be missing, which will guarantee 
#that these features are always in the output file, even with an 'NA' value. 
#'cutoff' applies a minimum read cutoff to the data. Any site below that cutoff is not counted.
#'site_cutoff_only', when set to True, will only apply the cutoff to the summing up the number of
#sites, but will calculated the total weighted methylation using all the data regardless of the cutoff.
def metaplot(mc_bed, gff, genome_file, output, mc_type=['CG','CHG','CHH'], flank_distance=2000, 
	window_number=20, primary_feature='gene', secondary_feature='CDS', chrs=[], feature_ids=[], 
	cutoff=0, site_cutoff_only=False):
	#List temporary files created to remove later
	files = ['temp.db', 'sf_bed.tmp', 'pf_bed.tmp', 'u_bed.tmp', 'd_bed.tmp', 'w_bed.tmp', 
		'p_bed.tmp', 'n_bed.tmp']
	#Use sf2pf_bed to map secondary_feature to primary_feature and create
	sf_pbt_bed,pf_list,pf_pbt_bed = sf2pf_bed(gff, primary_feature=primary_feature, 
		secondary_feature=secondary_feature, chrs=chrs, return_pf_bed=True)
	#If feature_ids are provided, then all analyses will be restricted to that list of ids
	if feature_ids:
		#Use 'name_filter' function to filter pf_pbt_bed and sf_pbt_bed for those features
		pf_pbt_bed = pf_pbt_bed.filter(name_filter,feature_ids).saveas('pf_bed_filtered.tmp')
		sf_pbt_bed = sf_pbt_bed.filter(name_filter,feature_ids).saveas('sf_bed_filtered.tmp')
		#Add to temporary file list
		files = files + ['sf_bed_filtered.tmp', 'pf_bed_filtered.tmp']
	#If a flanking distance is provided, create bedfiles for the upstream 'u_bed' and downstream 
	#'d_bed' regions
	if flank_distance:
		#Upstream
		u_pbt_bed = pbt.bedtool.BedTool.flank(pf_pbt_bed, g=genome_file, 
			l=flank_distance, r=0, s=True).saveas('u_bed.tmp')
		#Downstream
		d_pbt_bed = pbt.bedtool.BedTool.flank(pf_pbt_bed, g=genome_file, 
			l=0, r=flank_distance, s=True).saveas('d_bed.tmp')
		#Set the list of regions to include upstream, feature, & downstream
		regions=[u_pbt_bed,pf_pbt_bed,d_pbt_bed]
	else:
		#Otherwise the regions only correspond to the feature
		regions=[pf_pbt_bed]
	#Set windows to start at 1
	window_offset=0
	for f in regions:
		#Filter bed files based on strand. We do this, so that we can reverse the order for featuers 
		#on the negative strand, otherwise the order will be messed up
		#Positive strand features
		p_bed = f.filter(strand_filter,strand='+').saveas('p_bed.tmp')
		#Negative strand features
		n_bed = f.filter(strand_filter,strand='-').saveas('n_bed.tmp')
		#Make windows using pybedtools window_maker
		#Positive strand features
		pw_bed = pbt.bedtool.BedTool.window_maker(p_bed, b=p_bed, n=window_number, i='srcwinnum')
		#Now the negative strand features with the option 'reverse=True', this makes sures to put 
		#them in the correct order
		nw_bed = pbt.bedtool.BedTool.window_maker(n_bed, b=n_bed, n=window_number, i='srcwinnum', 
			reverse=True)
		#Combine windows and read in as a pandas dataframe so that we can modify the gene names
		wdf = pd.read_csv(pw_bed.cat(nw_bed, postmerge=False).fn, header=None, sep="\t")
		#Modify the fourth column with gene names, leaving only the window number
		wdf[3].replace('^.*_', '', inplace=True, regex=True)
		#Add the window_offset to the windows, that way it increases with each round
		wdf[3] = (wdf[3].astype(int) + window_offset).astype(str)
		#convert back to a pbt bedfile	
		#If we are on the primary feature bed, we intersect with the sf_pbt_bed to remove regions
		#not corresponding to the secondary feature. 
		if f == pf_pbt_bed:
			#Check if this is the first round, e.g. window_offset is equal to 0, if so we set it to
			#w_bed to equal it.
			if window_offset == 0:
				w_bed = pbt.bedtool.BedTool.intersect(pbt.BedTool.from_dataframe(wdf), 
					sf_pbt_bed).saveas('w_bed.tmp')
			#Otherwise, we concatinate it.
			else:
				w_bed.cat(pbt.bedtool.BedTool.intersect(pbt.BedTool.from_dataframe(wdf), 
					sf_pbt_bed), postmerge=False).saveas('w_bed.tmp')
		#If not pf_pbt_bed, then its the u_pbt_bed or d_pbt_bed. 
		else:
			#Check if this is the first round, e.g. window_offset is equal to 0, if so we set it to
			#w_bed to equal it.
			if window_offset == 0:
				w_bed = pbt.BedTool.from_dataframe(wdf).saveas('w_bed.tmp')
			#Otherwise, we concatinate it.
			else:
				w_bed.cat(pbt.BedTool.from_dataframe(wdf), postmerge=False).saveas('w_bed.tmp')
		#Add the window_number to the window_offset, which will then be used in the next iteration
		window_offset = window_offset + window_number
	#We now invoke map_to_bed to map out methylation data to our windows and write to output
	map_to_bed(mc_bed, w_bed, output, feature_ids=list(map(str, range(1, window_offset + 1))), 
		mc_type=mc_type, cutoff=cutoff, site_cutoff_only=site_cutoff_only)
	#And we delete our temporary gff database and bed file.
	for tmp_file in files:
		remove(tmp_file)


