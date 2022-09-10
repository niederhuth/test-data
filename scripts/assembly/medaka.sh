#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=800GB
#SBATCH --job-name medaka
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
input_dir="racon_*" #Directory or directories with fasta files to polish, assumes racon polished fasta!
model="r941_prom_high_g4011" #dna_r9.4.1_450bps_hac.cfg on guppy 4.2.2, guppy 4 models same as guppy 3.6
threads=20
datatype="ont"
batch_size=100

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
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)

#Extract reads from assembly job report
reads="$(grep reads: ${path2}/job_reports/${assembly}-*.SLURMout | head -1 | cut -d ' ' -f2)"

#Run medaka
for i in ${input_dir}
do
	path3="medaka_${i}"
	fasta="${i}/${i}.fa"
	if [ -d ${path3} ]
	then
		echo "Directory ${path3} exists."
		echo "To rerun medaka, delete ${path3} and resubmit."
	else
		echo "Running medaka on ${fasta}"
		medaka_consensus \
			-i ${path2}/${reads} \
			-d ${fasta} \
			-o ${path3} \
			-t ${threads} \
			-m ${model} \
			-b ${batch_size}
	fi
done

echo "Done"


