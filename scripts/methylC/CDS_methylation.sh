#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name CDS_methylation
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
datatype="methylC"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/methylC/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/methylC/lib:${LD_LIBRARY_PATH}"
#Export temporary directory paths
export TMPDIR=${PBS_O_WORKDIR}
export TMP=${PBS_O_WORKDIR}
export TEMP=${PBS_O_WORKDIR}

#Everything below this should not need changed, unless you are modifying the pipeline

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts\\/methylC\\/py/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path3=${datatype}

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $7}' \
	${path1}/samples.csv)

#Create output directory
if [ ! -d methylC/results ]
then
	mkdir methylC/results
fi

#Lets analyze the data
for i in ${genomes}
do
	#Set variables
	path4=$(pwd | sed s/${species}\\/.*/${species}\\/${i}/)
	version=$(ls ${path4}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	ref_path="${path4}/ref"
	output="methylC/results/${sample}_ref_${i}_CDS_methylation.tsv"
	#Analyze data
	if [ -f ${output} ]
	then
		echo "Existing allc file found, exiting"
		echo "To rerun this step, delete ${output} and resubmit"
	else
		echo "Calling methylation levels for gene coding sequences"
		python ${path2}/CDS_methylation.py \
			${sample} \
			${i} \
			${version} \
			${ref_path}
	fi
done

echo "Done"
