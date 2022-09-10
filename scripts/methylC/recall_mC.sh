#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name recall_mC
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
datatype="methylC"
threads=20
sort_mem="5G" #Memory used in sorting files
upstream=0 #Number of bases upstream, almost always 0
downstream=2 #Number of bases downstream, typically for plants 2, non-plants 1
SNPs=False #Add snp info? I really haven't had luck with this, so default is False
min_cov=3 #Minimum coverage necessary to call a base as methylated or not by binomial test
min_mapq=30 #Minimum mapping quality to keep read, default is 30

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/methylC/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/methylC/lib:${LD_LIBRARY_PATH}"

#Everything below this should not need changed, unless you are modifying the pipeline

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2=${datatype}

#Check if paired-end or single-end
path3="fastq/${datatype}"
if ls ${path3}/*_R1_001.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_R2_001.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		PE="True"
	else
		echo "Data is Single-end"
		PE="False"
	fi
elif ls ${path3}/*_1.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_2.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		PE="True"
	else
		echo "Data is Single-end"
		PE="False"
	fi
else
	echo "Data Missing"
fi

#Get list of genomes
read -a genomes <<< $(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $7}' \
	${path1}/samples.csv)

#Get unmethylated reference sequence. This is typically ChrL or ChrC
read -a unmethylated_control  <<<  $(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $10}' \
	${path1}/samples.csv)

#Recall methylation
for i in ${!genomes[@]}
do
	#Set variables
	path4=$(pwd | sed s/${species}\\/.*/${species}\\/${genomes[i]}/)
	version=$(ls ${path4}/ref/${genomes[i]}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	fasta="${path4}/ref/methylpy/${genomes[i]}-v${version}.fa"
	bam_file="${path2}/${sample}_ref_${genomes[i]}_processed_reads_no_clonal.bam"
	#Report unmethylated control
	echo " "
	echo "###############################"
	echo "Recalling methylation for alignment to ${genomes[i]}"
	echo "Unmethylated control chromosome is ${unmethylated_control[i]}"

	#Align Data
	if [ -f ${path2}/allc_${sample}_ref_${genomes[i]}.tsv.gz ]
	then
		echo "Existing allc file found, exiting"
		echo "To rerun this step, delete ${path2}/allc_${sample}_ref_${genomes[i]}.tsv.gz and resubmit"
	else
		methylpy call-methylation-state \
			--input-file ${bam_file} \
			--sample ${sample}_ref_${genomes[i]} \
			--ref-fasta ${fasta} \
			--paired-end ${PE} \
			--path-to-output ${path2} \
			--num-procs ${threads} \
			--num-upstream-bases ${upstream} \
			--num-downstream-bases ${downstream} \
			--remove-chr-prefix False \
			--add-snp-info ${SNPs} \
			--unmethylated-control ${unmethylated_control[i]} \
			--binom-test True \
			--min-mapq ${min_mapq} \
			--min-cov ${min_cov} 
	fi
done

echo "Done"




