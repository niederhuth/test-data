#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name methylpy
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
datatype="methylC"
threads=20
sort_mem="5G" #Memory used in sorting files
upstream=0 #Number of bases upstream, almost always 0
downstream=2 #Number of bases downstream, typically for plants 2, non-plants 1
aligner="bowtie2" #Can be bowtie2, bowtie, minimap2 -> Stongly recommend bowtie2
aligner_options="" #Specify aligner options. If empty, my own defaults will be set later in the script.
SNPs=False #Add snp info? I really haven't had luck with this, so default is Faalse
min_cov=3 #Minimum coverage necessary to call a base as methylated or not by binomial test
min_mapq=30 #Minimum mapping quality to keep read, default is 30
methylpy_trim_reads=True #Trim reads using cutadapt in methylpy (True) or before methylpy (False)
trimmer="cutadapt" #What tool to use for trimming (cutadapt or trimmomatic)
pbat=False #Were libraries prepared using post-bisulfite adapter tagging? default False

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/methylC/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/methylC/lib:${LD_LIBRARY_PATH}"
#Path to picard jar. Change this to appropriate path
picard="${conda}/envs/methylC/share/picard-2.25.1-0/"

#Everything below this should not need changed, unless you are modifying the pipeline

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2=${datatype}

#Get unmethylated reference sequence. This is typically chrL or chrC
unmethylated_control=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $10}' \
	${path1}/samples.csv)
echo "Unmethylated control is ${unmethylated_control}"

#Fastq files, these should not have to be changed, but should set automatically
path3="fastq/${datatype}"
r1="${path3}/combined.1.fastq.gz"
r2="${path3}/combined.2.fastq.gz"
t1="${path3}/trimmed.1.fastq.gz"
t2="${path3}/trimmed.2.fastq.gz"
t3="${path3}/trimmed.1.single.fastq.gz"
t4="${path3}/trimmed.2.single.fastq.gz"
if ls ${path3}/*_R1_001.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_R2_001.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_R1_001.fastq.gz > $r1
			cat ${path3}/*_R2_001.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_R1_001.fastq.gz > $r1
		fi
	fi
elif ls ${path3}/*_1.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_2.fastq.gz >/dev/null 2>&1	
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_1.fastq.gz > $r1
			cat ${path3}/*_2.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_1.fastq.gz > $r1
		fi
	fi
else
	echo "Data Missing"
fi

#Get read length
length=$(expr $(zcat ${r1} | head -2 | tail -1 | wc -c) - 1)

#Get adapters from samples.csv
adapters=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $8}' \
	${path1}/samples.csv)
#Path to adapter fastas 
adapter_path="${path1}/${adapters}_adapters.fa"

#Now get the adapters if using cutadapt or methylpy
#Need to keep this line regardless
if [ ${PE} = "TRUE" ]
then
	adapter1=$(grep -A 1 \> ${adapter_path} | grep -A 1 "1$" | tail -1)
	adapter2=$(grep -A 1 \> ${adapter_path} | grep -A 1 "2$" | tail -1)
else
	adapter_SE=$(grep -A 1 \> ${adapter_path} | tail -1)
fi

#If trimming adapters prior to methylpy
if [ ${methylpy_trim_reads} = "False" ]
then
	#Trim & QC reads
	if [ -f ${t1} ]
	then
		if [ ${PE} = "TRUE" ]
		then
			echo "To rerun this step, please delete ${t1} & ${t2} and resubmit"
		else
			echo "To rerun this step, please delete ${t1} and resubmit"
		fi
	else
		if [ ${PE} = "TRUE" ]
		then
			if [ ${trimmer} = "cutadapt" ]
			then
			echo "Running cutadapt"
			cutadapt \
				-j ${threads} \
				-a ${adapter1} \
				-A ${adapter2} \
				-q 10 \
				--trim-n \
				-m 30 \
				-e 0.2 \
				-O 3 \
				-n 1 \
				-o ${t1} \
				-p ${t2} \
				${r1} ${r2}
			elif [ ${trimmer} = "trimmomatic" ]
			then
				echo "Running trimmomatic"
				trimmomatic PE \
					-threads ${threads} \
					-phred33 \
					-trimlog ${path3}/trim_log.txt \
					-summary ${path3}/trim_summary.txt \
					${r1} ${r2} ${t1} ${t3} ${t2} ${t4} \
					ILLUMINACLIP:${adapter_path}:2:30:10:4:TRUE \
					MINLEN:30
			fi
			#If used Swift Biosciences Accel-NGS, then need to trim 15 bases
			if [ ${adapters} = "swift_mC-PE" ]
			then
				tmp1="${path3}/trimmed.1.fastq.gz"
				tmp2="${path3}/trimmed.2.fastq.gz"
				t1="${path3}/swift.1.fastq.gz"
				t2="${path3}/swift.2.fastq.gz"
				echo "Trimming adaptase added sequences"
				cutadapt \
					-j ${threads} \
					-u -15 \
					-u 4 \
					-U 17 \
					-l 100 \
					-m 30 \
					-o ${t1} \
					-p ${t2} \
					${tmp1} ${tmp2}
			fi
			echo "Running fastqc"
			mkdir ${path3}/fastqc
			fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${t2} ${r1} ${r2}
		elif [ ${PE} = "FALSE" ]
		then
			if [ ${trimmer} = "cutadapt" ]
			then
				echo "Running cutadapt"
				cutadapt \
					-j ${threads} \
					-a ${adapter_SE} \
					--trim-n \
					-m 30 \
					-e 0.1 \
					-O 3 \
					-n 1 \
					-o ${t1} \
					${r1}
			elif [ ${trimmer} = "trimmomatic" ]
			then
				echo "Running trimmomatic"
				trimmomatic SE \
					-threads ${threads} \
					-phred33 \
					-trimlog ${path3}/trim_log.txt \
					-summary ${path3}/trim_summary.txt \
					${r1} ${t1} \
					ILLUMINACLIP:${adapter_path}:2:30:10:4:TRUE \
					MINLEN:30 
			fi
			#If used Swift Biosciences Accel-NGS, then need to trim 15 bases
			if [ ${adapters} = "swift_mC-SE" ]
			then
				tmp1="${path3}/trimmed.1.fastq.gz"
				t1="${path3}/swift.1.fastq.gz"
				echo "Trimming adaptase added sequences"
				cutadapt \
					-j ${threads} \
					-u -15 \
					-u 4 \
					-m 30 \
					-o ${t1} \
					${tmp1}
			fi
			echo "Running fastqc"
			mkdir ${path3}/fastqc
			fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${r1}
		fi
	fi
else
	if [ ${PE} = "TRUE" ]
	then
		t1=${r1}
		t2=${r2}
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${r1} ${r2}
	else
		t1=${r1}
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${r1}
	fi
fi	

#Set aligner options
if [ -z ${aligner_options} ]
then
	if [ ${aligner} = "bowtie2" ]
	then
		aligner_options="--very-sensitive -X 1000"
	elif [ ${aligner} = "bowtie" ]
	then
		aligner_options=""
	elif [ ${aligner} = "minimap2" ]
	then
		aligner_options=""
	fi
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

#Create output directory
mkdir ${path2}

#Lets run methylpy
for i in ${!genomes[@]}
do
	#Set variables
	path4=$(pwd | sed s/${species}\\/.*/${species}\\/${genomes[i]}/)
	version=$(ls ${path4}/ref/${genomes[i]}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	f_ref="${path4}/ref/methylpy/${genomes[i]}-v${version}_f"
	r_ref="${path4}/ref/methylpy/${genomes[i]}-v${version}_r"
	fasta="${path4}/ref/methylpy/${genomes[i]}-v${version}.fa"
	#Report unmethylated control
	echo " "
	echo "###############################"
	echo "Aligning to genome ${genomes[i]}"
	echo "Unmethylated control chromosome is ${unmethylated_control[i]}"

	#Align Data
	if [ -f ${path2}/allc_${sample}.tsv.gz ]
	then
		echo "Existing allc file found, exiting"
		echo "To rerun this step, delete ${path2}/allc_${sample}_ref_${i}.tsv.gz and resubmit"
	else
		if [ ${PE} = "TRUE" ]
		then
			methylpy paired-end-pipeline \
				--read1-files ${t1} \
				--read2-files ${t2} \
				--sample ${sample}_ref_${genomes[i]} \
				--forward-ref ${f_ref} \
				--reverse-ref ${r_ref} \
				--ref-fasta ${fasta} \
				--libraries "libA" \
				--path-to-output ${path2} \
				--pbat ${pbat} \
				--num-procs ${threads} \
				--sort-mem ${sort_mem} \
				--num-upstream-bases ${upstream} \
				--num-downstream-bases ${downstream} \
				--trim-reads ${methylpy_trim_reads} \
				--aligner ${aligner} \
				--aligner-options "${aligner_options}" \
				--merge-by-max-mapq True \
				--remove-clonal True \
				--path-to-picard ${picard} \
				--adapter-seq-read1 ${adapter1} \
				--adapter-seq-read2 ${adapter2} \
				--remove-chr-prefix False \
				--add-snp-info ${SNPs} \
				--unmethylated-control ${unmethylated_control[i]} \
				--binom-test True \
				--min-mapq ${min_mapq} \
				--min-cov ${min_cov} \
				--max-adapter-removal 1 \
				--overlap-length 3 \
				--error-rate 0.1
		elif [ ${PE} = "FALSE" ]
		then
			methylpy single-end-pipeline \
				--read-files ${t1} \
	        		--sample ${sample}_ref_${genomes[i]} \
				--forward-ref ${f_ref} \
				--reverse-ref ${r_ref} \
				--ref-fasta ${fasta} \
				--libraries "libA" \
				--path-to-output ${path2} \
				--pbat ${pbat} \
				--num-procs ${threads} \
				--sort-mem ${sort_mem} \
				--num-upstream-bases ${upstream} \
				--num-downstream-bases ${downstream} \
				--trim-reads ${methylpy_trim_reads} \
				--aligner ${aligner} \
				--aligner-options "${aligner_options}" \
				--merge-by-max-mapq True \
				--remove-clonal True \
				--path-to-picard ${picard} \
				--adapter-seq ${adapter_SE} \
				--remove-chr-prefix False \
				--add-snp-info ${SNPs} \
				--unmethylated-control ${unmethylated_control[i]} \
				--binom-test True \
				--min-mapq ${min_mapq} \
				--min-cov ${min_cov} \
				--max-adapter-removal 1 \
				--overlap-length 3 \
				--error-rate 0.1
		fi
	fi
done

echo "Done"




