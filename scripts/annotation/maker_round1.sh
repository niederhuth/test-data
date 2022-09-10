#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200GB
#SBATCH --job-name maker_round1
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
SR=TRUE #Use short read transcript asemblies?
LR=TRUE #Use long read transcript assemblies?
fasta= #input fasta, if left blank, will look for it in current directory
blast_threads=1 #Leave 1 for MPI

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"
#Export path to agusutus config files
#export ZOE="${conda}/envs/maker" #Need to check
export AUGUSTUS_CONFIG_PATH="${conda}/envs/maker/config/"
#export REPEATMASKER_LIB_DIR=
#export REPEATMASKER_MATRICES_DIR=

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path3="maker_round1"

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

#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#Copy over and reformat repeatmasker gff
if [ ! -s repeatmasker.gff ]
then
	echo "Getting repeatmasker gff file"
	#Repeatmasker gff file is in gff2 and has non-unique IDs, which maker hates, so we have to fix this
	${conda}/envs/maker/share/RepeatMasker/util/rmOutToGFF3.pl ../repeatmasker/${fasta}.out | \
	sed s/RepeatMasker/repeatmasker/ | \
	sed s/dispersed_repeat/match/ | \
	grep -v \# | \
	sed s/Target\=/Repeat\=/ | \
 	awk -v OFS="\t" '{a+=1}{print $1,$2,$3,$4,$5,".",$7,$8,"ID="a";"$9}' > repeatmasker.gff
else
	echo "repeatmasker gff found"
fi

#Get protein gff
if [ ! -s protein_alignments_maker_input.gff ]
then
	echo "Getting protein alignments"
	#Get list of proteins to sources to use
	protein_sources="$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		'{if ($1 == a && $2 == b && $3 == c) print $4}' \
		${path1}/annotation/annotation_sources.csv)"
	for i in ${protein_sources}
	do
		echo "Adding ${i} protein alginments"
		cat ../exonerate/${i}/${i} >> protein_alignments
	done
	#Convert to gff
	perl ${path2}/annotation/pl/reformat_exonerate_protein_gff.pl \
		--input_gff protein_alignments \
		--output_gff tmp.gff
	#Sort gff file
	gff3_sort -g tmp.gff -og protein_alignments.gff
	rm tmp.gff
	#Reformatfor maker
	cat protein_alignments.gff | \
	sed 's/mRNA/protein_match/g' | \
	sed 's/exon/match_part/g' | \
	sed s/protein2genome/protein_gff\:protein2genome/ > protein_alignments_maker_input.gff
else
	echo "Protein alignment gff found"
fi

#Get est gff
if [ ! -s est_maker_input.gff ]
then
	echo "Getting est alignments"
	#Get list of short-read sources
	if [ ${SR} = TRUE ]
	then
		SR_sources="$(awk -v FS="," \
			-v a=${species} \
			-v b=${genotype} \
			-v c=${sample} \
			'{if ($1 == a && $2 == b && $3 == c) print $5}' \
			${path1}/annotation/annotation_sources.csv)"
	fi
	#Get list of long-read sources
	if [ ${LR} = TRUE ]
	then
		LR_sources="$(awk -v FS="," \
			-v a=${species} \
			-v b=${genotype} \
			-v c=${sample} \
			'{if ($1 == a && $2 == b && $3 == c) print $6}' \
			${path1}/annotation/annotation_sources.csv)"
	fi
	for i in ${SR_sources} ${LR_sources}
	do
		echo "Adding ${i} est alginments"
		cat ../stringtie/${i}_maker_input.gff >> tmp.gff
	done
	#Sort the gff file
	gff3_sort -g tmp.gff -og est_maker_input.gff
	rm tmp.gff
else
	echo "est gff found"
fi

#Run maker
echo "Running Maker Round 1 on ${fasta/.fa*/}"
maker \
	-q \
	-genome ../${fasta} \
	-cpus ${blast_threads} \
	${path1}/annotation/${path3}/*

#Get gff & fasta files
gff3_merge -d ${fasta/.f*/}.maker.output/${fasta/.f*/}_master_datastore_index.log
fasta_merge -d ${fasta/.f*/}.maker.output/${fasta/.f*/}_master_datastore_index.log

echo "Done"

