#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH --job-name polca
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
rounds=4
input="" #Can set to empty and script will find fasta in directory submitted

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/masurca/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/masurca/lib:$LD_LIBRARY_PATH"
#Export paths to MaSuRCA stuff
export PATH="$(pwd | sed s/data.*/scripts/)/assembly/MaSuRCA-4.0.5/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd | sed s/data.*/scripts/)/assembly/MaSuRCA-4.0.5/lib:$LD_LIBRARY_PATH"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/masurca/share/trimmomatic/adapters"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)

#Adapter fasta, set automatically from misc/samples.csv
adapters=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $8}' \
	${path1}/samples.csv)

#Fastq files, these should not have to be changed, but should set automatically
path3="${path2}/fastq/${datatype}"
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
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_R1_001.fastq.gz > $r1
			cat ${path3}/*_R2_001.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if [ -f ${t1} ]
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
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_1.fastq.gz > $r1
			cat ${path3}/*_2.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_1.fastq.gz > $r1
		fi
	fi
else
	echo "Data Missing"
fi

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
		echo "Running trimmomatic PE"
		trimmomatic PE \
			-threads ${threads} \
			-phred33 \
			-trimlog ${path3}/trim_log.txt \
			-summary ${path3}/trim_summary.txt \
			${r1} ${r2} ${t1} ${t3} ${t2} $t4 \
			ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:30
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${t2} ${r1} ${r2}
	elif [ ${PE} = "FALSE" ]
	then
		echo "Running trimmomatic SE"
		trimmomatic SE \
			-threads ${threads} \
			-phred33 \
			-trimlog ${path3}/trim_log.txt \
			-summary ${path3}/trim_summary.txt \
			${r1} ${t1} \
			ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:30 
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${r1}
	fi
fi
rm ${r1} ${r2}

#Define Read Group
ID=$(zcat ${t1} | head -1 | cut -d ':' -f 3,4 | tr ':' '.')
PU=$(zcat ${t1} | head -1 | cut -d ':' -f 3,4,10 | tr ':' '.')
SM=$(pwd | sed s/^.*\\///)
PL="ILLUMINA"
LB="lib1"

#Look for fasta file, there can only be one!
if [ -z ${input} ]
then
	echo "No input fasta provided, looking for fasta"
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
else
	echo "Input fasta: ${input}"
fi

#Loop over designated number of rounds for polishing
ref=../${input}
a=0
until [ ${a} -eq ${rounds} ]
do
	a=$(expr ${a} + 1)
	path4="polca_${a}"
	echo "Round ${a} of polishing" 
	if [ -d ${path4} ]							
	then
		echo "Directory ${path4} exists."
		echo "Checking files."
		cd ${path4}
	else
		mkdir ${path4}
		cd ${path4}
	fi
	#Polish with Polca
	if [ -s polca_${a}.fasta ]
	then
		echo "Round ${a} polishing found"
		echo "To rerun this step, delete ${path4}/polca_${a}.fasta and resubmit"
	else
		echo "Polishing data with Polca"
		if [ ${PE} = "TRUE" ]
		then
			polca.sh \
				-a ${ref} \
				-r '../../../fastq/wgs/trimmed.1.fastq.gz ../../../fastq/wgs/trimmed.2.fastq.gz' \
				-t ${threads}
		else
			polca.sh \
				-a ${ref} \
				-r '../../../fastq/wgs/trimmed.1.fastq.gz' \
				-t ${threads}
		fi
	fi
	mv ${input}.PolcaCorrected.fa polca_${a}.fasta
	ref="../${path4}/polca_${a}.fasta"
	input="polca_${a}.fasta"
	cd ../
	echo "Round ${a} of polishing complete"
done

echo "Done"

