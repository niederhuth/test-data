#!/bin/bash

OPTS=`getopt -a --longoptions chrUN,prefix:,justify:,zeros_at_end:,input_gff:,input_protein_fa:,input_transcript_fa:,output_prefix: -n "$0" -- "$@"`

eval set -- "${OPTS}"
while :
do
  case "$1" in
	--chrUN)				CHRUN=1			; shift   ;;
	--prefix)				PREFIX="$2"		; shift 2 ;;
	--justify)				JUSTIFY="$2"	; shift 2 ;;
	--zeros_at_end)			ZEROS="$2"		; shift 2 ;;
	--input_gff)			GFF="$2"		; shift 2 ;;
	--input_protein_fa)		PROTEINS="$2"	; shift 2 ;;
	--input_transcript_fa)	TRANSCRIPTS="$2"; shift 2 ;;
	--output_prefix)		OUTPUT="$2"		; shift 2 ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called.
    *) echo "Unexpected option: $1."
       usage ;;
	esac
done

##Arguments
#prefix=$1
#justify=$2
#zeros_at_end=$3
#input_gff=$4
#PROTEINS=$5
#TRANSCRIPTS=$6
#OUTPUT=$7

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"

#Get chromosome list
cut -f1 ${GFF} | grep -v \# | sort | uniq > chr_list

#Loop over chr_list and rename genes for each chromosome
mkdir tmp
cd tmp
cat ../chr_list | while read line
do
	echo "Working on ${line}"
	if [[ ${line:0:3} = "chr" ]]
	then
		awk -v a=${line} '$1==a' ../${GFF} > ${line}.gff
		if [ $(echo ${line} | wc -c) -gt 5 ]
		then
			chr=$(echo ${line} | sed s/chr//)
		else
			chr=$(echo ${line} | sed s/chr/0/)
		fi
		maker_map_ids \
		--prefix ${PREFIX}${chr}g \
		--justify ${JUSTIFY} \
		--iterate 1 \
		${line}.gff | awk -v a=${ZEROS} '{if ($2 ~ /-R/) print $0; else print $0a}' | sed s/-R/${ZEROS}./ > ${line}_renamed.map
	elif [[ CHRUN ]]
	then
		awk -v a=${line} '$1==a' ../${GFF} >> chrUN.gff
	else
		awk -v a=${line} '$1==a' ../${GFF} > ${line}.gff
		maker_map_ids \
			--prefix ${PREFIX}${line}g \
			--justify ${JUSTIFY} \
			--iterate 1 \
			${line}.gff | awk -v a=${ZEROS} '{if ($2 ~ /-R/) print $0; else print $0a}' | sed s/-R/${ZEROS}./ > ${line}_renamed.map
	fi
	chr=
done

if [[ chrUN ]]
then
	chr="UN"
	maker_map_ids \
		--prefix ${PREFIX}${chr}g \
		--justify ${JUSTIFY} \
		--iterate 1 \
		chrUN.gff | awk -v a=${ZEROS} '{if ($2 ~ /-R/) print $0; else print $0a}' | sed s/-R/${ZEROS}./ > chrUN_renamed.map
fi

#Combine files
cd ..
cat tmp/*_renamed.map > ${OUTPUT}-renamed-genes.map
rm -R tmp chr_list

#Rename gff & fasta files
map_gff_ids ${OUTPUT}-renamed-genes.map ${GFF}
map_fasta_ids ${OUTPUT}-renamed-genes.map ${TRANSCRIPTS}
map_fasta_ids ${OUTPUT}-renamed-genes.map ${PROTEINS}

echo "Done"
