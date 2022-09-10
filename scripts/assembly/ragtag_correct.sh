#!/bin/bash --login
#SBATCH --time=128:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name ragtag_correct
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
read_validation=TRUE #TRUE or FALSE
reads= #if read validation TRUE and this is blank, will look for reads based on datatype
threads=20
datatype="wgs"
aligner="unimap"
min_len=1000
merge_dist=100000
break_dist=5000
cov_win_size=8000
min_cov=
max_cov=

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

#Look for reads
if [ ${read_validation} = "TRUE" ]
then
	if [ ${datatype} = "ont" ]
	then
		read_type="ont"
		echo "${path2}/fastq/${datatype}/clean.fastq.gz" > reads.fofn
	elif [ ${datatype} = "wgs" ]
	then
		read_type="sr"
		if [ -f ${path2}/fastq/${datatype}/trimmed.2.fastq.gz ]
		then
			echo "${path2}/fastq/${datatype}/trimmed.1.fastq.gz" > reads.fofn
			echo "${path2}/fastq/${datatype}/trimmed.2.fastq.gz" >> reads.fofn
		elif [ -f ${path2}/fastq/${datatype}/trimmed.1.fastq.gz ]
		then
			echo "${path2}/fastq/${datatype}/trimmed.1.fastq.gz" > reads.fofn
		fi
	fi
fi

#Check for coverage limits
if [ -z ${min_cov} ]
then
	if [ -z ${max_cov} ]
	then
		cov_vars="-v ${cov_win_size}"
	else
		cov_vars="-v ${cov_win_size} --max-cov ${max_cov}"
	fi
else
	if [ -z ${max_cov} ]
	then
		cov_vars="-v ${cov_win_size} --min-cov ${min_cov}"
	else
		cov_vars="-v ${cov_win_size} --min-cov ${min_cov} --max-cov ${max_cov}"
	fi
fi

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
	${path1}/samples.csv)

#Run ragtag correct
for i in ${genomes}
do
	#Find reference genome
	species2=$(echo ${i} | sed s/_.*//)
	genotype2=$(echo ${i} | sed s/${species2}_// | sed s/_.*//)
	path3=$(pwd | sed s/data\\/.*/data/)
	path4="${path3}/${species2}/${genotype2}"
	version=$(ls ${path4}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	ref="${path4}/ref/${i}-v${version}.fa"
	echo "Running ragtag correct with ${i}-v${version}.fa as reference genome"
	
	#Check if reads are provided
	if [ ${read_validation} = "FALSE" ]
	then
		#If not, then run based on just alignment to the genome
		echo "No reads provided, running without validation"
		ragtag.py correct \
			-t ${threads} \
			-o ${i}_ragtag_correct \
			--aligner ${aligner} \
			-R ${min_len} \
			--remove-small \
			-d ${merge_dist} \
			-b ${break_dist} \
			-u \
			${ref} ${input}
	else
		#If reads provided, run with read validation
		echo "Using ${reads} for validation"
		ragtag.py correct \
			-t ${threads} \
			-o ${i}_ragtag_correct \
			--aligner ${aligner} \
			-R ${min_len} \
			--remove-small \
			-d ${merge_dist} \
			-b ${break_dist} \
			-u \
			--read-aligner minimap2 \
			-F reads.fofn \
			-T ${read_type} \
			${cov_vars} \
			${ref} ${input}
	fi
done

echo "Done"


