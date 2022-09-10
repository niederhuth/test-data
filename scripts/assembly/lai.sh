#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200GB
#SBATCH --job-name LAI
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
engine="rmblast" #crossmatch, wublast, abblast, ncbi, rmblast, hmmer
threads=10 #Actual cores used by RM are thread number multiplied by the cores for search engine used.
			#These are: RMBlast=4 cores, ABBlast=4 cores, nhmmer=2 cores, crossmatch=1 core
			#So 10 threads with RMBlast actually needs 40 cores!
#Set what repeat library to use. This is currently set to a set of denovo TEs identified by EDTA
LTR_lib="../edta/*.fa.mod.EDTA.raw/LTR/*.fa.mod.LTRlib.fa" #non-redundant LTRlib from LTR_retriever
LTR_list="../edta/*.fa.mod.EDTA.raw/LTR/*.fa.mod.pass.list"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/EDTA/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/EDTA/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path2="lai"

#Look for fasta file, there can only be one!
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

#mask with repeatmasker
if ls ${path2} >/dev/null 2>&1
then
	echo "Previous LAI results found."
	echo "Quitting. To repeat this analysis, delete ${path2} and resubmit."
else
	mkdir ${path2}
	cd ${path2}
	cp ${LTR_lib} ./LTRlib.fa
	cp ${LTR_list} ./LTR_list.pass.list
	echo "Running RepeatMasker"
	RepeatMasker \
		-e ${engine} \
		-pa ${threads} \
		-q \
		-norna \
		-no_is \
		-div 40 \
		-cutoff 225 \
		-dir .
		-lib LTRlib.fa \
		../${input}
	echo "Running LAI"
	LAI \
		-genome ../${input} \
		-intact LTR_list.pass.list \
		-all ${input}.out
fi

echo "Done"
