#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name busco
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
datatype="ont"
remove_tmp_files="yes" #if yes, delete a lot of the busco files to save space, if no, keep them
lineage="" #if left empty, will check misc/samples to find lineage

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/busco/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/busco/lib:$LD_LIBRARY_PATH"
#Export path to agusutus config files
export AUGUSTUS_CONFIG_PATH="${conda}/envs/busco/config"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2="busco"
output="${genotype}_${assembly}"

#Check if lineage is set, otherwise, get lineage to use for Busco analysis
if [ -z ${lineage} ]
then
	lineage=$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		-v d=${condition} \
		-v e=${datatype} \
		'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $10}' \
		${path1}/samples.csv)
fi
echo "Using lineage dataset: ${lineage}"

#Look for fasta file, there can only be one!
if ls *.fa >/dev/null 2>&1
then
	input="../"$(ls *fa | sed s/.*\ //)
	echo "Fasta file ${input} found"
elif ls *.fasta >/dev/null 2>&1
then
	input="../"$(ls *fasta | sed s/.*\ //)
	echo "Fasta file ${input} found"
elif ls *.fna >/dev/null 2>&1
then
	input="../"$(ls *fna | sed s/.*\ //)
	echo "Fasta file ${input} found"
else
	echo "No fasta file found, please check and restart"
fi

#Run busco
if [ -d ${path2} ]
then
	echo "Previous Busco run found, restarting from last startpoint"
	echo "To restart from scratch, delete the directory ${assembly}/busco/${output} and resubmit"
	cd ${path2}
	busco \
		-c ${threads} \
		-m genome \
		-i ${input} \
		-l ${lineage} \
		-o ${output} \
		--restart
else
	echo "Running Busco analysis for ${species} ${genotype} ${sample} ${assembly}"
	mkdir ${path2}
	cd ${path2}
	busco \
		-c ${threads} \
		-m genome \
		-i ${input} \
		-l ${lineage} \
		-o ${output}
fi

if [ ${remove_tmp_files} = "yes" ]
then
	rm -R */*/*output */*/busco_sequences/ busco_downloads/ */logs/
fi

echo "Done"
