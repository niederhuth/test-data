#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name medaka_stich
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#If medaka_consensus fails after producing consensus_probs.hdf, use this to finish pipeline

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
input="consensus_probs.hdf"
threads=20
output="consensus.fasta"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
dir=$(pwd | sed s/^.*\\///)

#Run medaka
fasta="../${dir/medaka_/}/${dir/medaka_/}.fa"
echo "Running medaka stitch for ${fasta}"
medaka stitch \
	${input} \
	${fasta} \
	${output} \
	--threads ${threads} \

echo "Done"

