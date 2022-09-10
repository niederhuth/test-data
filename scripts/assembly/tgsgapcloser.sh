#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name tgsgapcloser
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20 #doesn't seem to want to use more than 6
datatype="ont" #ont or pb
min_identity=0.3 #minimum identity for filter reads: 0.3 for ont by default/0.2 for pb by default.
min_match=300 #min match length for filter reads: 300bp for ont by default/200bp for pb by default.
chunk_number=10 #number of chunks to use
racon=TRUE
racon_rounds=1
pilon=FALSE
ngs="wgs/" #Have not implemented yet!
pilon_round=2
pilon_mem="300G" #memory used for pilon , 300G for default.
input= #input fasta, if left blank, will look for it in current directory, mutually exclusive with input_dir

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)
reads="${path2}/fastq/${datatype}/clean.fastq.gz"

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

#Set options
options="--thread ${threads} --tgstype ${datatype}"
if [ ${racon} = "TRUE" ]
then
	options="${options} --racon ./ --r_round ${racon_rounds} --chunk ${chunks}"
fi
if [ ${pilon} = "TRUE" ]
then
	options="${options} --ngs ${ngs} --pilon ./ --samtools ./ --java ./ --p_round ${pilon_round} --pilon_mem ${pilon_mem} --chunk ${chunks}"
fi
if [[ ${racon} = "FALSE" && ${pilon} = "FALSE" ]]
then
	options="--ne"
fi

#Make and cd to workign directory
if [ -d tgsgapcloser ]
then
	cd tgsgapcloser
else
	mkdir tgsgapcloser
	cd tgsgapcloser
fi

#Run tgsgapcloser
echo "Running tgsgapcloser"
tgsgapcloser \
	--min_idy ${min_identity} \
	--min_match ${min_match} \
	--scaff ../${input} \
	--reads ${reads} \
	--output tgsgapcloser \
	${options}

echo "Done"

