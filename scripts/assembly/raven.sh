#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name raven
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
iterations=2

#This should match the dataype in the misc/samples.csv file
datatype="ont" 

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
path2="raven"
reads="fastq/${datatype}/clean.fastq.gz"

#Declare reads
echo "reads: ${reads}"

#Get genome size estimate
genomeSize=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $9}' \
	${path1}/samples.csv)

#Run raven
if [ -d ${path2} ]
then
	echo "Previous raven assembly detected, restarting from previous run"
	echo "To start from the beginning, please delete the directory ${path2} and resubmit"
	cd ${path2}
	raven \
		--polishing-rounds ${iterations} \
		--threads ${threads} \
		--resume \
		../${reads} > assembly.fasta
else
	mkdir ${path2}
	cd ${path2}
	echo "Beginning raven assembly"
	raven \
		--polishing-rounds ${iterations} \
		--threads ${threads} \
		../${reads} > assembly.fasta 
fi

echo "Done"
