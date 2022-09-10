#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name decontaminate
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"
database="gtdb-rs202.genomic.k31.sbt.json"

#Set variables
#Currently none to set

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
microbes="$(pwd | sed s/data.*/data/)/microbe-database/${database}"
path2="contaminants"

#Check for microbe-genbank/genbank-k31.sbt.json
if [ -s ${microbes} ]
then
	echo "genbank-k31.sbt.json  found"
else
	echo "Could not find genbank-k31.sbt.json"
	echo "Run scripts/assembly/download_microbe_database.sh and resubmit"
fi

#Look for fasta file, there can only be one!
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

#Download database
mkdir ${path2}
echo "Checking for contaminants with gather-by-contig.py"
${path1}/py/gather-by-contig.py \
	${input} \
	${microbes} \
	--output-nomatch ${path2}/clean.fa \
	--output-match ${path2}/putative_contaminants.txt \
	--csv ${path2}/summary.csv

echo "Done"
