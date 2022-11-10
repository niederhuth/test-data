#!/bin/bash --login
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name CDS-GC-content
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta= #Input fasta of seqs, if left blank will search for it in the reference for that genotype
datatype="cds-primary" #cds/cds-primary, cds = all transcripts, cds-primary = primary transcripts only

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/orthofinder/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/orthofinder/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/scripts/)
path2=$(pwd | sed s/data.*/data/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition=""
path3="GC_content"

#Make and cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Find input CDS sequence
if [ -z ${fasta} ]
then
	version=$(ls ${path2}/${species}/${genotype}/ref/annotations/${genotype}-v${version}-${datatype}.fa | \
		sed s/.*\-v// | sed s/\-${datatype}.fa)
	fasta="${path2}/${species}/${genotype}/ref/annotations/${genotype}-v${version}-${datatype}.fa"
fi

#Set output
output=${genotype}-v${version}-${datatype}_GC_content.tsv

#get total weighted mC
echo "Get GC content of ${genotype}-v${version}-${datatype}.fa"
python ${path1}/sequence/py/CDS_GC_content.py ${fasta} ${output}

echo "Done"