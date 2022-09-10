#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name purge_haplotigs_step2
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
low=4
mid=16
high=60
suspect=80
junk=60
align_cov=70
max_match=250
fasta=""

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
assembly=$(pwd | sed s/^.*\\///)
path1=$(pwd | sed s/${genotype}.*/${genotype}/)
path2="purge_haplotigs"

#Output location
echo "Purging Duplicates for ${species} ${genotype} ${sample} ${assembly}"

#Look for fasta file, there can only be one!
if [ -z ${fasta} ]
then
	echo "No input fasta provided, looking for fasta"
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
else
	echo "Input fasta: ${fasta}"
fi

if [ -d ${path2} ]
then
	cd ${path2}
else
	echo "${path2} directory not found. Has purge_haplotigs_step1 been run?"
	echo "If not, first run purge_haplotigs_step1.sh and then rerun this step."
fi

#Get coverage stats
if [ -s coverage_stats.csv ]
then
	echo "Coverage stats found, proceeding to haplotig purging."
	echo "To repeat this step, delete ${path2}/coverage_stats.csv and resubmit."
else
	echo "Generating coverage statistics"
	purge_haplotigs cov \
		-in aligned.bam.gencov \
		-low ${low} \
		-high ${high} \
		-mid ${mid} \
		-out coverage_stats.csv \
		-junk ${junk} \
		-suspect ${suspect}
fi

#Purge duplicates
if [ -s curated.fasta ]
then
	echo "Aligned reads found, proceeding to coverage statistics."
	echo "To repeat this step, delete ${path2}/curated.fasta and resubmit."
else
	echo "Purging haplotigs"
	purge_haplotigs purge \
		-genome ../${fasta} \
		-coverage coverage_stats.csv \
		-threads ${threads} \
		-outprefix curated \
		-dotplots \
		-bam aligned.bam \
		-align_cov ${align_cov} \
		-max_match ${max_match} 
fi

echo "Done"
