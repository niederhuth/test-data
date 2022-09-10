#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200GB
#SBATCH --job-name train_augustus
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
maker_dir="../maker_round1" #maker dir with gff & transcripts, not needed if input_gff & transcripts specified
input_gff= #input gff file, if left blank will look in maker_dir
transcripts= #transcript fasta file, if left blank will look in maker_dir
AUGUSTUS_SPECIES_NAME="" #Species name for Augustus, if left blank will automatically set
fasta= #input fasta, if left blank, will look for it in current directory

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"
#Export path to agusutus config files
export AUGUSTUS_CONFIG_PATH="${conda}/envs/maker/config/"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path3="augustus_training"

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

#Set some more inputs
if [ -z ${AUGUSTUS_SPECIES_NAME} ]
then
	AUGUSTUS_SPECIES_NAME="${species}-${genotype}-${sample}"
fi
if [ -z ${input_gff} ]
then
	input_gff="${maker_dir}/${fasta/.f*/}.all.gff"
fi
if [ -z ${transcripts} ]
then
	transcripts="${maker_dir}/${fasta/.f*/}.all.maker.transcripts.fasta"
fi


#Run maker2zff
#For determining which genes are High Confidence for Retraining, there are 7 criteria.
#-c fraction  The fraction of splice sites confirmed by an EST alignment, default 0.5
#-e fraction  The fraction of exons that overlap an EST alignment, default 0.5
#-o fraction  The fraction of exons that overlap any evidence (EST or Protein), default 0.5
#-a fraction  The fraction of splice sites confirmed by an ab-initio prediction, default 0
#-t fraction  The fraction of exons that overlap an ab-initio prediction, default 0
#-l number    The min length of the protein sequence produced by the mRNA
#-x number    Max AED to allow 0.5 is default
#-n           No filtering.  Accept all.
echo "Running maker2zff"
maker2zff \
	-c 0.5 \
	-e 0.5 \
	-o 0.5 \
	-a 0 \
	-t 0 \
	-l 200 \
	-x 0.2 \
	${input_gff}

#Run fathom
#-validate [-quiet]
#-gene-stats [-errors-ok -nucleotide -dinucleotide]
#-categorize <int>
#-export <int> [-plus -errors-ok]
#-extract <feature> -length <int> -offset <int>
#-exon-intron
#-split <-number <int> | -training <float> | -GC <float> | -repeat <float>>
#-ace-format <-gene-method <string> [-dna -extra <string>]>
#-compare-genes <predictions> [-details]
#-score-genes <hmm> [-errors-ok]
#-filter-genes <hmm> -min-score <float> -min-length <int>
echo "Running fathom"
fathom \
	genome.ann \
	genome.dna \
	-categorize 1000

NUMFOUND="`grep -c '>' uni.ann`"

if [ ${NUMFOUND} -gt 599 ]
then
	NUMFOUND=600
fi

TEMPSPLIT=$((NUMFOUND/2))
NUMSPLIT=${TEMPSPLIT/.*}

echo "number found after fathom: ${NUMFOUND}"
echo "number after split: ${NUMSPLIT}"

#Convert the uni.ann and uni.dna output from fathom into a genbank formatted file.
#fathom_to_genbank.pl is from https://github.com/Childs-Lab/GC_specific_MAKER.
echo "Conveting fathom to genbank format"
perl ${path2}/annotation/pl/fathom_to_genbank.pl \
	--annotation_file uni.ann \
	--dna_file uni.dna \
	--genbank_file augustus.gb \
	--number ${NUMFOUND}

#Get the subset of fastas that correspond to the genes in the genbank file.
#get_subset_of_fastas.pl is from github: https://github.com/Childs-Lab/GC_specific_MAKER.
echo "Making genbank gene list"
perl -e  'while (my $line = <>){ if ($line =~ /^LOCUS\s+(\S+)/) { print "$1\n"; } }' \
augustus.gb > genbank_gene_list.txt

echo "Subsetting fastas"
perl ${path2}/annotation/pl/get_subset_of_fastas.pl \
	-l genbank_gene_list.txt \
	-f uni.dna \
	-o genbank_gene_seqs.fasta

#Split genes into test and training files.
echo "Randomly splitting genes into test and training files"
randomSplit.pl augustus.gb ${NUMSPLIT}

#We will use autoAug.pl for training because we have transcript alignments that will be used as hints.
#The etraining and optimize_augustus.pl scripts from AUGUSTUS will not be used as they do not allow 
#for the use of hints.
#rm -R /mnt/home/niederhu/miniconda3/envs/maker/config//species/Mguttatus-S1-S1
echo "Running first round of autoAug.pl"
autoAug.pl \
	--species=${AUGUSTUS_SPECIES_NAME} \
	--genome=genbank_gene_seqs.fasta \
	--trainingset=augustus.gb.train \
	--cdna=${transcripts} \
	--noutr

#Run the batch scripts generated by autoAug.pl
echo "Running first round of batch scripts"
cd autoAug/autoAugPred_abinitio/shells
x=1
while [ -e ./aug${x} ]
do
    echo "A.  $x"
    ./aug${x}
    let x=x+1
done

#Run the next command as indicated by autoAug.pl in step 5.
#When above jobs are finished, continue by running the command autoAug.pl
cd ../../../
echo "Running second round of autoAug.pl"
autoAug.pl \
	--species=${AUGUSTUS_SPECIES_NAME} \
	--genome=genbank_gene_seqs.fasta \
	--useexisting \
	--hints=autoAug/hints/hints.E.gff \
	-v -v -v --index=1

#Run the batch scripts generated by autoAug.pl
echo "Running second round of batch scripts"
cd autoAug/autoAugPred_hints/shells/
let x=1
while [ -e ./aug${x} ]
do
    echo "B.  $x"
    ./aug${x}
    let x=x+1
done

#Checked sensitivity and specificity of the newly trained AUGUSTUS HMM by using the test data.
cd ../../../
echo "Checking sensitivity and specificity"
augustus \
	--species=${AUGUSTUS_SPECIES_NAME} \
	augustus.gb.test

echo "Done"

