#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100GB
#SBATCH --job-name filter_genes
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
hmmscan_evalue='1e-5' #Cutoff evalue for hmmscan against pfam
blast_evalue='1e-10' #Cutoff evalue for blast
fasta= #genome fasta, if left blank will search for in submission directory
gff= #input gff, if left blank will search in maker_dir
transcripts= #input transcripts fa, if left blank will search in maker_dir
proteins= #input transcripts fa, if left blank will search in maker_dir
maker_dir="../maker_round2"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path3="gene_filtering"

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

#Set various datasets
if [ -z ${gff} ]
then
	gff=${maker_dir}/${fasta/.fa/}.all.gff
fi
if [ -z ${transcripts} ]
then
	transcripts=${maker_dir}/${fasta/.fa/}.all.maker.transcripts.fasta
fi
if [ -z ${proteins} ]
then
	proteins=${maker_dir}/${fasta/.fa/}.all.maker.proteins.fasta
fi

#Download Pfam A hmm database
echo "Downloading Pfam-A"
wget -q http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

#Prepare Pfam A hmm-database
echo "Preparing Pfam A hmm database"
hmmpress Pfam-A.hmm

#Search Pfam A hmm domains
echo "Searching against Pfam A hmm profiles"
hmmscan \
	--domE ${hmmscan_evalue} \
	-E ${hmmscan_evalue} \
	--cpu ${threads} \
	-o pfam_alignments.out \
	--tblout prot_domains.out \
	Pfam-A.hmm \
	${proteins}

#Prefilter genes for those without any support (AED=1, or non-sig pfam result)
#Generate maker standard gene list
echo "Generating prefiltered list of mRNAs"
perl ${path2}/annotation/pl/generate_maker_standard_gene_list.pl \
	--input_gff ${gff} \
	--pfam_results prot_domains.out \
    --pfam_cutoff 1e-10 \
	--output_file ${fasta/.fa/}_prefilter_mRNA.txt \
	> bad_mRNAs.txt

echo "Making prefiltered gff file"
perl ${path2}/annotation/pl/create_maker_standard_gff.pl \
	--maker_standard_gene_list ${fasta/.fa/}_prefilter_mRNA.txt \
	--input_gff ${gff} \
	--output_gff ${fasta/.fa/}.gff

echo "Making prefiltered transcripts fasta"
perl ${path2}/annotation/pl/get_subset_of_fastas.pl \
	-l ${fasta/.fa/}_prefilter_mRNA.txt \
	-f ${transcripts} \
	-o ${fasta/.fa/}-transcripts.fa
	
echo "Making prefiltered protein fasta"
perl ${path2}/annotation/pl/get_subset_of_fastas.pl \
	-l ${fasta/.fa/}_prefilter_mRNA.txt \
	-f ${proteins} \
	-o ${fasta/.fa/}-proteins.fa

#Download and make Transposase blast DB
echo "Downloading Tpases020812"
wget -q http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz
gunzip Tpases020812.gz
#remove trailing white space in the Tpases020812 file
sed -i 's/[ \t]*$//' Tpases020812
#Make diamond blast DB
echo "Making Transposase diamond blast DB"
diamond makedb \
	--threads ${threads} \
	--in Tpases020812 \
	--db Tpases020812.dmnd

#Run diamond blastp against Transposases 
echo "Running diamond blastp on Transposases"
diamond blastp \
	--threads ${threads} \
	--db Tpases020812.dmnd \
	--query ${proteins} \
	--out TE_blast.out\
	--evalue ${blast_evalue} \
	--outfmt 6

#Eliminated gypsy hmmscan as this resulted in large amount of false positives
#Download Gypsy DB hmm files and format the hmm database
#echo "Downloading GyDB_collection"
#wget https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
#unzip GyDB_collection.zip
#echo "Combining GyDB hmm profiles"
#cat GyDB_collection/profiles/*hmm > all_gypsy.hmm
#echo "Formatting all_gypsy.hmm database"
#hmmpress all_gypsy.hmm

#Search gypsy hmm profiles
#echo "Searching against gypsy hmm profiles"
#hmmscan \
#	--domE 1e-5 \
#	-E 1e-5 \
#	--cpu ${threads} \
#	-o gypsy_alignments.out \
#	--tblout gypsyHMM_analysis.out \
#	all_gypsy.hmm \
#	${proteins}

#Create a genelist with no TEs
echo "Creating gene list with TEs removed"
python ${path2}/annotation/py/create_TE_noTE_genelist.py \
	--input_geneList ${fasta/.fa/}_prefilter_mRNA.txt \
	--pfamhmm prot_domains.out \
	--TEpfam_list ${path1}/annotation/TE_Pfam_domains.txt \
	--TEblast TE_blast.out \
	--output_non_TE_genes noTE_mRNA_list.txt \
	--output_TE_like_genes potential_TE_like_mRNA_list.txt
	#--TEhmm gypsyHMM_analysis.out #Eliminated gypsy hmmscan as this resulted in large amount of false positives

#Filter out potential TEs
echo "Making no TE gff"
perl ${path2}/annotation/pl/create_maker_standard_gff.pl \
	--maker_standard_gene_list noTE_mRNA_list.txt \
	--input_gff ${fasta/.fa/}.gff \
	--output_gff ${fasta/.fa/}-noTE.gff

echo "Making no TE Transcripts fasta"
perl ${path2}/annotation/pl/get_subset_of_fastas.pl \
	-l noTE_mRNA_list.txt \
	-f ${fasta/.fa/}-transcripts.fa \
	-o ${fasta/.fa/}-transcripts-noTE.fa

echo "Making no TE proteins fasta"
perl ${path2}/annotation/pl/get_subset_of_fastas.pl \
	-l noTE_mRNA_list.txt \
	-f ${fasta/.fa/}-proteins.fa \
	-o ${fasta/.fa/}-proteins-noTE.fa

#Convert mRNA id to gene id in noTE and potential_TE_like lists
sed s/-mRNA.*// noTE_mRNA_list.txt > noTE_gene_list.txt
sed s/-mRNA.*// potential_TE_like_mRNA_list.txt > potential_TE_like_gene_list.txt

echo "Done"

