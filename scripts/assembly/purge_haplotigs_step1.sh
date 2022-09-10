#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200GB
#SBATCH --job-name purge_haplotigs_step1
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
threads2=20
datatype="ont"
depth=200
fasta=""

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
assembly=$(pwd | sed s/.*${genotype}\\/${sample}\\/// | sed s/\\/.*//)
path1=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)
path2="purge_haplotigs"

#Output location
echo "Purging Duplicates for ${species} ${genotype} ${sample} ${assembly}"

#Extract reads from assembly job report
reads="$(grep reads: ${path1}/job_reports/${assembly}-*.SLURMout | head -1 | cut -d ' ' -f2)"

#Change preset based on datatype
if [ ${datatype} = "ont" ]
then
	preset="map-ont"
elif [ ${datatype} = "ont-cor" ]
then
	preset="map-ont"
elif [ ${datatype} = "pac" ]
then
	preset="map-pb"
elif [ ${datatype} = "pac-cor" ]
then
	preset="map-pb"
elif [ ${datatype} = "hifi" ]
then
	preset="map-pb"
else
	echo "Do not recognize ${datatype}"
	echo "Please check and resubmit"
fi

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

#Create output directory and change directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Align reads to assembly
if [ -s aligned.sam ]
then
	echo "Aligned reads found, proceeding to bam conversion & sorting."
	echo "To repeat this step, delete ${path2}/aligned.sam and resubmit."
else
	echo "Aligning reads to assembly"
	minimap2 \
		-a \
		-t ${threads} \
		-x ${preset} \
		../${fasta} \
		${path1}/${reads} > aligned.sam
fi

if [ -s aligned.bam ]
then
	echo "Bam file found, proceeding to read-depth histogram."
	echo "To repeat this step, delete ${path2}/aligned.bam and resubmit."
else
	echo "Sorting and converting to bam"
	samtools view -@ ${threads2} -bSh aligned.sam | samtools sort -@ ${threads2} > aligned.bam
	echo "Indexing aligned.bam"
	samtools index aligned.bam
fi

#Generate read-depth histogram
if [ -s aligned.bam.genecov ]
then
	echo "Read depth histogram found."
	echo "To repeat this step, delete ${path2}/aligned.bam.genecov and resubmit."
else
	echo "Generating read-depth histogram"
	purge_haplotigs hist \
		-bam aligned.bam \
		-genome ../${fasta} \
		-threads ${threads} \
		-depth ${depth}
fi

echo "purge_haplotigs_step1 complete. Choose read-depth cutoffs and proceed to purge_haplotigs_step_2"
echo "Done"

