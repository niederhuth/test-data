#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=100GB
#SBATCH --job-name align-rna-SR-annotation
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40 #for HISAT2
threads2=4 #for samtools sort
options="--rna-strandness R --very-sensitive" #R/RF for dUTP fr-firststrand F/FR for fr-secondstrand
fasta= #input fasta, if left blank, will look for it in current directory

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/transcript-assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/transcript-assembly/lib:${LD_LIBRARY_PATH}"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/transcript-assembly/share/trimmomatic/adapters"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="annotation"
datatype="rna"
assembly=$(pwd | sed s/^.*\\///)
path2="SRrna"

#Look for fasta file, there can only be one!
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

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Make HISAT2 index
if [ ! -d hisat2_index ]
then
	mkdir hisat2_index
fi

if [ ! -s hisat2_index/${fasta}.8.ht2 ]
then
	echo "Making HISAT2 index"
	hisat2-build ../${fasta} hisat2_index/${fasta}
else
	echo "HISAT2 index found"
fi

#Get list of datasets
datasets=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	'{if ($1 == a && $2 == b && $3 == c) print $5}' \
	${path1}/annotation/annotation_sources.csv)

#Iterate over each dataset, trim, and align
for i in ${datasets}
do
	species2=$(echo ${i} | sed s/_.*//)
	genotype2=$(echo ${i} | sed s/${species2}_// | sed s/_.*//)
	sample2=$(echo ${i} | sed s/.*_//)
	path3=$(pwd | sed s/data\\/.*/data/)
	path4="${path3}/${species2}/${genotype2}/${sample2}/fastq/${datatype}"
	bam="${species2}_${genotype2}_${sample2}.bam"

	echo "Working on data ${species2} ${genotype2} ${sample2}"

	#Adapter fasta, set automatically from misc/samples.csv
	adapters=$(awk -v FS="," \
		-v a=${species2} \
		-v b=${genotype2} \
		-v c=${sample2} \
		-v d=${condition} \
		-v e=${datatype} \
		'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $8}' \
		${path1}/samples.csv)

	#Fastq files, these should not have to be changed, but should set automatically
	r1="${path4}/combined.1.fastq.gz"
	r2="${path4}/combined.2.fastq.gz"
	t1="${path4}/trimmed.1.fastq.gz"
	t2="${path4}/trimmed.2.fastq.gz"
	t3="${path4}/trimmed.1.single.fastq.gz"
	t4="${path4}/trimmed.2.single.fastq.gz"
	if ls ${path4}/*_R1_001.fastq.gz >/dev/null 2>&1
	then
		if ls ${path4}/*_R2_001.fastq.gz >/dev/null 2>&1
		then
			echo "Data is Paired-end"
			PE="TRUE"
			if ls ${t1} >/dev/null 2>&1
			then
				echo "Trimmed reads found, skipping trimming"
			else
				cat ${path4}/*_R1_001.fastq.gz > $r1
				cat ${path4}/*_R2_001.fastq.gz > $r2
			fi
		else
			echo "Data is Single-end"
			echo "Single-end ${datatype}? If this is wrong, double check and restart"
			PE="FALSE"
			if ls ${t1} >/dev/null 2>&1
			then
				echo "Trimmed reads found, skipping trimming"
			else
				cat ${path4}/*_R1_001.fastq.gz > $r1
			fi
		fi
	elif ls ${path4}/*_1.fastq.gz >/dev/null 2>&1
	then
		if ls ${path4}/*_2.fastq.gz >/dev/null 2>&1	
		then
			echo "Data is Paired-end"
			PE="TRUE"
			if ls ${t1} >/dev/null 2>&1
			then
				echo "Trimmed reads found, skipping trimming"
			else
				cat ${path4}/*_1.fastq.gz > $r1
				cat ${path4}/*_2.fastq.gz > $r2
			fi
		else
			echo "Data is Single-end"
			echo "Single-end ${datatype}? If this is wrong, double check and restart"
			PE="FALSE"
			if ls ${t1} >/dev/null 2>&1
			then
				echo "Trimmed reads found, skipping trimming"
			else
				cat ${path4}/*_1.fastq.gz > $r1
			fi
		fi
	elif ls ${path4}/SRR*.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Single-end"
		echo "Single-end ${datatype}? If this is wrong, double check and restart"
		PE="FALSE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path4}/SRR*.fastq.gz > $r1
		fi
	else
		echo "Data Missing"
	fi

	#Trim & QC reads
	if ls ${t1} >/dev/null 2>&1
	then
		if [ ${PE} = "TRUE" ]
		then
			echo "To rerun this step, please delete ${t1} & ${t2} and resubmit"
			fastq="${t1} ${t2}"
		else
			echo "To rerun this step, please delete ${t1} and resubmit"
			fastq=${t1}
		fi
	else
		if [ ${PE} = "TRUE" ]
		then
			echo "Running trimmomatic PE"
			trimmomatic PE \
				-threads ${threads} \
				-phred33 \
				-trimlog ${path4}/trim_log.txt \
				-summary ${path4}/trim_summary.txt \
				${r1} ${r2} ${t1} ${t3} ${t2} ${t4} \
				ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
				MINLEN:30
			echo "Running fastqc"
			mkdir ${path4}/fastqc
			fastqc -t ${threads} -o ${path4}/fastqc/ ${t1} ${t2} ${r1} ${r2}
			#Note that I am currently ignoring any unpaired reads
			fastq="${t1} ${t2}"
		elif [ ${PE} = "FALSE" ]
		then
			echo "Running trimmomatic SE"
			trimmomatic SE \
				-threads ${threads} \
				-phred33 \
				-trimlog ${path4}/trim_log.txt \
				-summary ${path4}/trim_summary.txt \
				${r1} ${t1} \
				ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
				MINLEN:30 
			echo "Running fastqc"
			mkdir ${path4}/fastqc
			fastqc -t ${threads} -o ${path4}/fastqc/ ${t1} ${r1}
			fastq=${t1}
		fi
	fi

	#Align with HISAT2
	if [ -s ${bam} ]
	then
		echo "${bam} found. Skipping alignment."
		echo "To rerun this step, please delete ${bam} and resubmit."
	else
		if [ ${PE} = "TRUE" ]
		then
			echo "Running HISAT2 for ${species2} ${genotype2} ${sample2} against ${fasta}"
			hisat2 ${options} \
				-p ${threads} \
				--dta \
				--no-mixed \
				--no-discordant \
				-x hisat2_index/${fasta} \
				-1 ${t1} \
				-2 ${t2} | \
				samtools view -@ ${threads2} -bSh | \
				samtools sort -@ ${threads2} > ${bam}
		elif [ ${PE} = "FALSE" ]
		then
			echo "Running HISAT2 for ${species2} ${genotype2} ${sample2} against ${fasta}"
			hisat2 ${options} \
				-p ${threads} \
				--dta \
				-x hisat2_index/${fasta} \
				-U ${t1} | \
				samtools view -@ ${threads2} -bSh | \
				samtools sort -@ ${threads2} > ${bam}		
		fi
	fi
done

echo "Done"
