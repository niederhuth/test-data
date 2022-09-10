#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200GB
#SBATCH --job-name alignment-stats
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)

#Fastq files, these should not have to be changed, but should set automatically
path3="${path2}/fastq/${datatype}"
t1="${path3}/trimmed.1.fastq.gz"
t2="${path3}/trimmed.2.fastq.gz"

#Look for fasta file, there can only be one!
if [ -z ${input} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		input=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		input=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fna >/dev/null 2>&1
	then
		input=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${input} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${input}"
fi

#check for index, if absent, make it
if [ -s ${input}.sa ]
then
	echo "BWA index found"
else
	echo "Indexing fasta"
	bwa index ${input}
fi

#Align reads
echo "Aligning reads to ${input} with bwa mem"
bwa mem -t ${threads} -M ${input} ${t1} ${t2} | samtools view -@ 4 -bSh | samtools sort -@ 4 > wgs.bam

#Index bam file
echo "Indexing bam file"
samtools index wgs.bam

#Get flagstats
echo "Getting flagstats"
samtools flagstat wgs.bam > alignment.flagstat

#Cleanup
rm wgs.bam wgs.bam.bai

echo "Done"

