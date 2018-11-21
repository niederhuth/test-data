import pybedtools as pbt
import pandas as pd
import numpy as np
import itertools

#interpret sequence context, taken from methylpy utilities
def expand_nucleotide_code(mc_type=["C"]):
	iub_dict = {"N":["A","C","G","T","N"],
	"H":["A","C","T","H"],
	"D":["A","G","T","D"],
	"B":["C","G","T","B"],
	"A":["A","C","G","A"],
	"R":["A","G","R"],
	"Y":["C","T","Y"],
	"K":["G","T","K"],
	"M":["A","C","M"],
	"S":["G","C","S"],
	"W":["A","T","W"],
	"C":["C"],
	"G":["G"],
	"T":["T"],
	"A":["A"]}
	mc_class = list(mc_type) # copy
	if "C" in mc_type:
		mc_class.extend(["CGN", "CHG", "CHH","CNN"])
	elif "CG" in mc_type:
		mc_class.extend(["CGN"])
	mc_class_final = []
	for motif in mc_class:
		mc_class_final.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in motif])])
	return(set(mc_class_final))

#function for filtering annotation files based on feature (gene, exon, mRNA, etc)
def annotation_filter(x,feature):
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
	#read in allc file to pandas dataframe
	a = pd.read_table(allc,dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int})
	#add new columns
	a['pos2'] = a.pos
	a['name'] = a.index
	a['score'] = "."
	#reorder columns
	a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
	#return bed file
	return a

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

#create filtered allc of sites mapping to annotations
def allc_annotation_filter(allc,annotations,genome_file,updown_stream=2000,first_feature=(),second_feature=(),filter_chr=[]):
	#
    bed = pbt.BedTool(annotations).filter(feat_filter,first_feature).filter(chr_filter,filter_chr)
	#
    flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=updown_stream,r=updown_stream,s=True).saveas('f_tmp')
	#
    cds_bed = pbt.BedTool(features).filter(feat_filter,second_feature).filter(chr_filter,filter_chr).saveas('c_tmp')
	#
    bed = cds_bed.cat(flank_bed, postmerge=False)
	#
    mC_bed = allc2bed(allc)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
	#
    m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
    m = m.drop_duplicates()
    m.to_csv('CDS_allc.tmp', sep='\t', index=False)
    os.remove('f_tmp')
    os.remove('c_tmp')

#
def metaplot(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],window_number=1,updown_stream=2000,cutoff=0,filter_chr=[]):
	a = allc2bed(allc)
	#create output data frame
	c = []
	columns=["","'Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC'"]
	for d in mc_type:
		for e in columns:
			c = c + [d + "_" + e]
	b = pd.DataFrame(columns=["Window"]+c)

	#list regions to iterate over: U:upstream, B=feature body, D=downstream
	regions=["U","B","D"]
	#if updownstream set to 0, only region to look at is B=feature body
	if updown_stream = 0:
		regions=["B"]
	#set window number to 1
	window = 1
	#iterate over each region and collect methylation data
	for i in regions:
		#filter bed files based on strand
		p_bed = bed.filter(strand_filter,strand='+').saveas('p_tmp')
		n_bed = bed.filter(strand_filter,strand='-').saveas('n_tmp')
		#make windows
		pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=window_number,i='srcwinnum')
		nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,n=window_number,i='srcwinnum',reverse=True)
		#combine windows
		w_bed = pw_bed.cat(nw_bed, postmerge=False)
		#intersect bedfiles with pybedtools
		allc_mapping = pbt.bedtool.BedTool.intersect(mC_bed,intersect_bed,wa=True,wb=True)
		#convert to pandas dataframe
		m = pd.read_table(allc_mapping.fn,usecols=[13,6,7,8,9])
		#split srcwinnum
		h = m[13].str.split("_", n = 1, expand = True)
		#make new columns from srcwinnum
		m["Name"] = h[0]
		m["Window"] = h[1]
		#drop column 13
		m.drop(columns = [13], inplace = True)[["Name","Window",6,7,8,9]]
		#iterate list of window numbers
		for x in range(1,window_number+1):
			#count windows
			window = window + 1
			#filter for rows matching specific window number
			c = a[a.chr.isin(x)]
			#make list for methylation data
			e = [window]
			#iterate over each mC type and run get_mC_data
			for z in mc_type:
				get_mC_data(c,mc_type=Z,cutoff=cutoff)
				#Calculate weighted methylation
				g = g + [np.float64(d[4])/np.float64(d[3])]
				#add this to list of data for other mc_types
				e = e + g
				#append the results for that window to the dataframe
				b = b.append(pd.DataFrame([e],columns=columns),ignore_index=True)
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b

#plot methylation levels for genes
def gene_metaplot(allc,features,genome_file,output=(),ignoreStrand=False,windows=60,
                  updown_stream=2000,cutoff=0,first_feature=(),second_feature=(),
                  filter_chr=[],remove_tmp=True):
    map2features(allc,features,genome_file,updown_stream,first_feature,second_feature,filter_chr)
    feature_metaplot('CDS_allc.tmp',features,genome_file,output,ignoreStrand,
                     windows,updown_stream,cutoff,first_feature,filter_chr)
    if remove_tmp:
        os.remove('CDS_allc.tmp')
