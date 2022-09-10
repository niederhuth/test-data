#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=220GB
#SBATCH --job-name kmers
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=10
mer_length=21
hash_size="20G"
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
path2="fastq/${datatype}/jellyfish"

#Create output directory
mkdir ${path2}
cd ${path2}

#Count kmers
echo "Counting kmers with jellyfish"
jellyfish count \
	-m ${mer_length} \
	-s ${hash_size} \
	-t ${threads} \
	-C <(zcat ../trimmed.1.fastq.gz) <(zcat ../trimmed.2.fastq.gz)

#Generate histogram
echo "Creating kmer histogram"
jellyfish histo \
	-t ${threads} \
	 mer_counts.jf > reads.histo

#Run Genomescope 2.0
echo "Running GenomeScope 2.0"
genomescope2 \
	-i reads.histo \
	-o output_dir \
	-k ${mer_length}

echo "Done"
