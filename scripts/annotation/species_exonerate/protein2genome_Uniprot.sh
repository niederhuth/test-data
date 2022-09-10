#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name protein2genome-Uniprot
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta=$(ls -l repeatmasker/*.fa.masked | sed s/.*\ //)
proteins="proteins/Uniprot-plant.fa"
query_chunks=5
target_chunks=10
minintron=10
maxintron=5000
bestn=5
ryo=">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path2="exonerate"

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Loop over protein sources and run exonerate
for i in ${proteins}
do
	echo "Working on ${i} proteins"
	outdir=$(echo ${i} | sed s/.*\\/// | sed s/\.fa//)
	a=1
	b=1
	if [ -d ${outdir} ]
	then
		echo "${outdir} directory present, checking files"
	else
		mkdir ${outdir}
	fi
	while [ ${a} -le ${target_chunks} ]
	do
		while [ ${b} -le ${query_chunks} ]
		do
			if [ -s ${outdir}/target_chunk_${a}_query_chunk_${b} ]
			then
				if [ $(wc -l ${outdir}/target_chunk_${a}_query_chunk_${b} | cut -d ' ' -f1) -gt 20 ]
				then
					echo "${outdir} target_chunk_${a}_query_chunk_${b} already complete" 
					echo "Skipping to next chunk"
					b=$(expr ${b} + 1)
				else
					rm ${outdir}/target_chunk_${a}_query_chunk_${b}
				fi
			fi
			if [ ! -s ${outdir}/target_chunk_${a}_query_chunk_${b} ]
			then
				if [ ${b} -le ${query_chunks} ]
				then
					echo "Aligning ${outdir} target_chunk_${a}_query_chunk_${b} on ${fasta}"
					exonerate \
						--model protein2genome \
						--bestn ${bestn} \
						--minintron ${minintron} \
						--maxintron ${maxintron} \
						--querychunkid ${b} \
						--querychunktotal ${query_chunks} \
						--targetchunkid ${a} \
						--targetchunktotal ${target_chunks} \
						--query ../${i} \
						--target ../${fasta} \
						--showtargetgff yes \
						--showalignment no \
						--showvulgar no \
						--ryo "${ryo}" > ${outdir}/target_chunk_${a}_query_chunk_${b}
					b=$(expr ${b} + 1)
				fi
			fi
		done
		b=1
		a=$(expr ${a} + 1)
	done
done

echo "Done"
