#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name stringtie
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
#To mimic --conservative set min_multi_exon_reads=1.5, min_iso_frac=0.05, trim=FALSE
threads=20
reformat_for_maker=TRUE #Convert gtf to gff for maker
SRread=TRUE #Run stringtie on short-read data
LRread=TRUE #Run stringtie on long-read data
SR_read_type="rf" #fr: fr-secondstrand, rf: fr-firststrand (dUTP method)
trim=TRUE #use coverage based trimming of transcript ends
min_multi_exon_reads="1" #min reads per bp cov for multi-exon transcript default: 1
min_single_exon_reads_SR="4.75" #min reads per bp cov for single-exon transcript from SR data, default: 4.75 for short reads
min_single_exon_reads_LR="1.5" #min reads per bp cov for single-exon transcript from LR data, default: 1.5 for long-reads
min_iso_frac="0.01" #minimum isoform fraction default: 0.01
max_gap_SR="50" #maximum gap allowed between read mappings default for SR data: 50 for short reads
max_gap_LR="0" #maximum gap allowed between read mappings default for LR data: 0 for long-reads
min_transcript_len="200" #minimum assembled transcript length default: 200
min_anchor_len="10" #minimum anchor length for junctions default: 10
min_junc_cov="1" #minimum junction coverage default: 1
frac_multi_hit="1" #fraction of bundle allowed to be covered by multi-hit reads default: 1
LR_splice_window="25" #window around possibly erroneous splice sites from long reads default: 25
combine_bams=FALSE #use all bam files for same run or separate runs, right now I do not recommend doing this

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/transcript-assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/transcript-assembly/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="annotation"
assembly=$(pwd | sed s/^.*\\///)
path2="stringtie"

#Make output dir and cd
if [ ! -d ${path2} ]
then
	mkdir ${path2}
fi
cd ${path2}

#Get list of datasets
if [ ${SRread} = TRUE ]
then
	SR_datasets="$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		'{if ($1 == a && $2 == b && $3 == c) print $5}' \
		${path1}/annotation/annotation_sources.csv)"
	#Loop over and assemble list of bam files
	for i in ${SR_datasets}
	do
		SR_bam_list="${SR_bam_list} ../SRrna/${i}.bam"
	done
	echo ${SR_bam_list}
fi
if [ ${LRread} = TRUE ]
then
	LR_datasets="$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		'{if ($1 == a && $2 == b && $3 == c) print $6}' \
		${path1}/annotation/annotation_sources.csv)"
	#Loop over and assemble list of bam files
	for i in ${LR_datasets}
	do
		LR_bam_list="${LR_bam_list} ../LRrna/${i}.bam"
	done
	echo ${LR_bam_list}
fi

#Add various settings
settings="-v -p ${threads}"
if [ ${trim} = FALSE ]
then
	settings="${settings} -t"
fi
#Run stringtie with bam files combined
if [ ${combine_bams} = TRUE ]
then
	echo "combine_bams is turned on"
	if [[ ${SRread} = TRUE && ${LRread} = TRUE ]]
	then
		echo "Using both short & long reads"
		settings="${SR_bam_list} ${LR_bam_list} ${settings} --mix -o combined.gtf"
	fi
	if [[ ${SRread} = TRUE && ${LRread} = FALSE ]]
	then
		echo "Using only short reads"
		if [ ${SR_read_type} = "rf" ]
		then
			settings="${SR_bam_list} ${settings} --rf  -g ${max_gap_SR} -s ${min_single_exon_reads_SR} -o SR_combined.gtf" 
		elif [ ${SR_read_type} = "fr" ]
		then
			settings="${SR_bam_list} ${settings} --fr  -g ${max_gap_SR} -s ${min_single_exon_reads_SR} -o SR_combined.gtf"
		fi
	fi
	if [[ ${SRread} = FALSE && ${LRread} = TRUE ]]
	then
		echo "Using only long reads"
		settings="${LR_bam_list} ${settings} -L -E ${LR_splice_window} -g ${max_gap_LR} -s ${min_single_exon_reads_LR} -o LR_combined.gtf"
	fi
	#Run stringtie
	echo ${settings}
	echo "Running stringtie on combined bam files"
	stringtie \
		${settings} \
		-c ${min_multi_exon_reads} \
		-f ${min_iso_frac} \
		-m ${min_transcript_len} \
		-a ${min_anchor_len} \
		-j ${min_junc_cov} \
		-M ${frac_multi_hit} \
		-l ${output}
fi
#Run stringtie on each bam file separately
if [ ${combine_bams} = FALSE ]
then
	echo "combine_bams is turned off"
	if [ ${SRread} = TRUE ]
	then
		echo "Running stringtie on short reads"
		if [ ${SR_read_type} = "rf" ]
		then
			settings="${settings} --rf"
		elif [ ${SR_read_type} = "fr" ]
		then
			settings="${settings} --fr"
		fi
		for i in ${SR_bam_list}
		do
			echo ${settings}
			output=$(echo ${i} | sed s/.*SRrna\\/// | sed s/.bam//)
			#Run stringtie
			echo "Running Stringtie on ${output}"
			stringtie ${i} \
			${settings} \
			-c ${min_multi_exon_reads} \
			-s ${min_single_exon_reads_SR} \
			-f ${min_iso_frac} \
			-g ${max_gap_SR} \
			-m ${min_transcript_len} \
			-a ${min_anchor_len} \
			-j ${min_junc_cov} \
			-M ${frac_multi_hit} \
			-l ${output} \
			-o ${output}.gtf
		done
	fi
	if [ ${LRread} = TRUE ]
	then
		echo "Running stringtie on long reads"
		settings="${settings} -L -E ${LR_splice_window}"
		for i in ${LR_bam_list}
		do
			output=$(echo ${i} | sed s/.*LRrna\\/// | sed s/.bam//)
			#Run stringtie
			echo "Running Stringtie on ${output}"
			stringtie ${i} \
			${settings} \
			-c ${min_multi_exon_reads} \
			-s ${min_single_exon_reads_LR} \
			-f ${min_iso_frac} \
			-g ${max_gap_LR} \
			-m ${min_transcript_len} \
			-a ${min_anchor_len} \
			-j ${min_junc_cov} \
			-M ${frac_multi_hit} \
			-l ${output} \
			-o ${output}.gtf
		done
	fi
fi

#Convert to gff3 and reformat for maker
if [ ${reformat_for_maker} = TRUE ]
then
	for i in *gtf
	do
		output2=$(echo ${i} | sed s/.gtf//)
		#Convert gtf to gff3
		echo "Converting ${i} to gff3"
		gffread ${i} -o tmp.gff

		#Sort the gff file. Maybe unnecessary, but just in case
		echo "Sorting gff file"
		gff3_sort -g tmp.gff -og ${output2}.gff
		rm tmp.gff

		#Modify gff for maker
		echo "Modifying for maker"
		cat ${output2}.gff | sed 's/transcript/expressed_sequence_match/g' | \
		sed 's/exon/match_part/g' | sed s/StringTie/est_gff\:est2genome/ > ${output2}_maker_input.gff
	done
fi

echo "Done"

