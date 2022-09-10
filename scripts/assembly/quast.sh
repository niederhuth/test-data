#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name quast
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
ref="IM62/ref/IM62-v2.fa"
features="IM62/ref/annotations/IM62-v2.gff"
eukaryote="true" #true or false, will modify arguments
fragmented="true" #true or false, will modify arguments
datatype="ont"
reads="fastq/ont/clean.fastq.gz"
ambiquity_score="" #between 0.8-1.0
min_contig=500

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/quast/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/quast/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
path1=$(pwd | sed s/${genotype}.*/${genotype}/)
path2=$(pwd | sed s/${species}.*/${species}/)

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

#Set arguments
args="--threads ${threads}" 

#Set arguments based on above variables
if [ ! -z ${min_contig} ]
then
	args="${args} --min-contig ${min_contig}"
fi

#Set arguments if eukaryote
if [ ${eukaryote} = "true" ]
then
	echo "Genome is eukaryotic"
	args="${args} --eukaryote --large --k-mer-stats"
fi

#Set arguments if fragmented
if [ ${fragmented} = "true" ]
then
	echo "Genome is fragmented"
	args="${args} --fragmented"
fi

#If reference genome is provided, add to args
if [ ! -z ${ref} ]
then
	echo "Reference genome provided"
	args="${args} -r ${path2}/${ref}"
fi

#If gene annotations from reference is provided, add to args
if [ ! -z ${features} ]
then
	echo "Reference annotations provdided"
	args="${args} -g ${path2}/${features}"
fi

#Change args based on datatype
if [ ${datatype} = "ont" ]
then
	if [ -z ${reads} ]
	then
		echo "No reads found"
                echo "If this is a mistake, check and resubmit"
	else
		args="${args} --nanopore ${path1}/${reads} --upper-bound-assembly --est-insert-size 255"
	fi
else
	echo "Do not recognize ${datatype}"
	echo "Please check and resubmit"
fi

#Run Quast
echo "Running Quast"
quast \
	${args} \
	${input}

echo "Done"

