import os
import sys
import re
import pandas as pd

functionsfile=re.sub('data.*','scripts/methylC/py/functions.py',os.getcwd())
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions as functions

#define variables
allc='methylC/allc_'+sys.argv[1]+'_ref_'+sys.argv[2]+'.tsv.gz'
gff=sys.argv[4]+'/annotations/'+sys.argv[2]+'-v'+sys.argv[3]+'.gff'
genome_file=sys.argv[4]+'/methylpy/'+sys.argv[2]+'-v'+sys.argv[3]+'.fa.fai'
filter_chr=['ChrL','ChrC','ChrM']
mc_type=['CG','CHG','CHH','CH']
flank_distance=2000
window_number=20
cutoff=3
site_cutoff_only=False
primary_feature='gene'
secondary_feature='CDS'
bedfile=sys.argv[1]+'_ref_'+sys.argv[2]+"_allc.bed"

#Check if a feature_ids are provided, if so load it as feature_ids
if sys.argv[5]:
	feature_ids = list(pd.read_csv(sys.argv[6], header=None, usecols=[0], dtype='str', sep="\t")[0])
	output='methylC/results/'+sys.argv[1]+'_'+sys.argv[5]+'_ref_'+sys.argv[2]+'_metaplot.tsv'
else:
	feature_ids=[]
	output='methylC/results/'+sys.argv[1]+'_ref_'+sys.argv[2]+'_metaplot.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file, header=None, usecols=[0], dtype='str', sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#Convert allc to bed file
print('Converting allc file to bed file')
functions.allc2bed(allc,bedfile)

#get gene methylation data
print('Getting gene methylation data')
functions.metaplot(bedfile, gff, genome_file, output, mc_type=mc_type, flank_distance=flank_distance, 
	window_number=window_number, primary_feature=primary_feature, secondary_feature=secondary_feature, 
	chrs=chrs, feature_ids=feature_ids, cutoff=cutoff, site_cutoff_only=site_cutoff_only)

#Delete the allc bedfile to save space.
print('Cleaning up')
os.remove(bedfile)

