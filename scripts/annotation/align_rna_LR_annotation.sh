#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=100GB
#SBATCH --job-name align-rna-LR-annotation
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40 #for minimap2
threads2=4 #for samtools sort
options="--secondary=no" #Additional options
input_data="*.fasta.gz" 
fasta= #input fasta, if left blank, will look for it in current directory

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/transcript-assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/transcript-assembly/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="annotation"
assembly=$(pwd | sed s/^.*\\///)
path2="LRrna"

#Look for fasta file, there can only be one!
if ls *.fa >/dev/null 2>&1
then
	fasta=$(ls *fa | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fasta >/dev/null 2>&1
then
	fasta=$(ls *fasta | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fna >/dev/null 2>&1
then
	fasta=$(ls *fna | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
else
	echo "No fasta file found, please check and restart"
fi

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Get list of datasets
datasets=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	'{if ($1 == a && $2 == b && $3 == c) print $6}' \
	${path1}/annotation/annotation_sources.csv)

#Iterate over each dataset, trim, and align
for i in ${datasets}
do
	species2=$(echo ${i} | sed s/_.*//)
	genotype2=$(echo ${i} | sed s/${species2}_// | sed s/_.*//)
	sample2=$(echo ${i} | sed s/${species2}_${genotype2}_// | sed s/_.*//)
	datatype=$(echo ${i} | sed s/.*_//) 
	path3=$(pwd | sed s/data\\/.*/data/)
	path4="${path3}/${species2}/${genotype2}/${sample2}/fastq/${datatype}"
	reads="${path4}/${input_data}"
	bam="${species2}_${genotype2}_${sample2}_${datatype}.bam"

	echo "Working on data ${species2} ${genotype2} ${sample2} ${datatype}"

	#Set various options
	options="-t ${threads} ${options}"
	if [ ${datatype} = "pb-rna" ]
	then
		options="${options} -ax splice:hq -uf"
		echo "Datatype is ${datatype}, setting options to minimap2 ${options}"
	elif [ ${datatype} = "ont-direct-rna" ]
	then
		options="${options} -ax splice -uf -k14"
		echo "Datatype is ${datatype}, setting options to minimap2 ${options}"
	elif [ ${datatype} = "ont-2D-rna" ]
	then
		options="${options} -ax splice"
		echo "Datatype is ${datatype}, setting options to minimap2 ${options}"
	else
		echo "Datatype not recognized"
	fi

	#Align with minimap2
	if [ -s ${bam} ]
	then
		echo "${bam} found. Skipping alignment."
		echo "To rerun this step, please delete ${bam} and resubmit."
	else
		echo "Aligning reads with minimap2"
		minimap2 \
			${options} \
			../${fasta} \
			${reads} | \
			samtools view -@ ${threads2} -bSh | \
			samtools sort -@ ${threads2} > ${bam}
	fi
done

echo "Done"
