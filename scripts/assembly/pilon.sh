#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH --job-name pilon
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
rounds=4
fix="all" #,breaks,novel"
java_options="-Xmx490G"
input="" #Can set to empty and script will find fasta in directory submitted
long_read=FALSE #TRUE or FALSE, from experience, illumina only tends to work better
long_read_type="ont"
threads2=5

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:$LD_LIBRARY_PATH"
#Path to picard
picard="${conda}/envs/polishing/share/picard-*/picard.jar"
#Path to pilon
pilon="${conda}/envs/polishing/share/pilon-*/pilon.jar"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/polishing/share/trimmomatic/adapters"

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
	path4="pilon_${a}"
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
	#Index fasta
	if [ -s ${ref}.sa ]
	then
		echo "BWA index found"
	else
		echo "Indexing fasta"
		bwa index ${ref}
	fi
	#Align data with bwa mem
	if [ -s round_${a}.bam ]
	then
		echo "Existing bam file found, skipping to mark duplicates"
		echo "To rerun this step, delete round_${a}.bam and resubmit"
	else
		echo "Running BWA for ${sample} to ${i} genome"
		if [ ${PE} = "TRUE" ]
		then
			bwa mem -t ${threads} \
				-R "@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}\tPU:${PU}" \
				-M ${ref} ${t1} ${t2} | \
				samtools view -@ 4 -bSh | samtools sort -@ 4 > round_${a}.bam
		elif [ ${PE} = "False" ]
		then	
			bwa mem -t ${threads} \
				-R "@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}\tPU:${PU}" \
				-M ${ref} ${t1} | samtools view -@ 4 -bSh | samtools sort -@ 4 > round_${a}.bam
		fi
		echo "Indexing round_${a}.bam"
		samtools index round_${a}.bam
	fi
	#Mark Duplicates
	if [ -s round_${a}_md.bam ]
	then
		echo "Existing marked duplicate bam found"
		echo "To rerun this step, delete round_${a}_md.bam and resubmit"
		bam_files="--frags round_${a}_md.bam"
	else	
		echo "Marking duplicates for round ${a}"
		java ${java_options} -jar ${picard} MarkDuplicates \
			-I round_${a}.bam \
			-O round_${a}_md.bam \
			-M round_${a}_md_metrics.txt
		echo "Indexing round ${a} marked duplicate bam"
		samtools index round_${a}_md.bam
		bam_files="--frags round_${a}_md.bam"
		#Alignment Stats
		echo "Getting round ${a} marked duplicate alignment stats"
		samtools flagstat round_${a}_md.bam > round_${a}_md.bam.flagstats
	fi
	if [ ${long_read} = "TRUE" ]
	then
		#Change preset based on datatype
		if [ ${long_read_type} = "ont" ]
		then
			preset="map-ont"
			bam_files="${bam_files} --nanopore long_reads.bam"
		elif [ ${long_read_type} = "ont-cor" ]
		then
			preset="map-ont"
			bam_files="${bam_files} --nanopore long_reads.bam"
		elif [ ${long_read_type} = "pac" ]
		then
			preset="map-pb"
			bam_files="${bam_files} --pacbio long_reads.bam"
		elif [ ${long_read_type} = "pac-cor" ]
		then
			preset="map-pb"
			bam_files="${bam_files} --pacbio long_reads.bam"
		elif [ ${long_read_type} = "hifi" ]
		then
			preset="map-pb"
			bam_files="${bam_files} --pacbio long_reads.bam"
		else
			echo "Do not recognize ${long_read_type}"
			echo "Please check and resubmit"
		fi
		if [ -s long_reads.bam ]
		then
			echo "Existing marked duplicate bam found"
			echo "To rerun this step, delete long_reads.bam and resubmit"
			bam_files="${bam_files} --nanopore long_reads.bam"
		else
			echo "Aligning ${long_read_type} reads to assembly"
			minimap2 \
				-a \
				-t ${threads} \
				-x ${preset} \
				${ref} \
				${path2}/fastq/${long_read_type}/clean.fastq.gz > long_reads.sam 
			echo "Sorting and converting to bam"
			samtools view -@ ${threads2} -bSh long_reads.sam | \
			samtools sort -@ ${threads2} > long_reads.bam
			echo "Indexing aligned.bam"
			samtools index long_reads.bam
			bam_files="${bam_files} --nanopore long_reads.bam"
		fi
	fi
	#Polish with Pilon
	if [ -s pilon_${a}.fasta ]
	then
		echo "Round ${a} polishing found"
		echo "To rerun this step, delete ${path4}/pilon_${a}.fasta and resubmit"
	else
		echo "Polishing data with Pilon"
		java ${java_options} -jar ${pilon} \
			--genome ${ref} \
			${bam_files} \
			--diploid \
			--fix ${fix} \
			--output pilon_${a} \
			--changes \
			--tracks
	fi
	ref="../${path4}/pilon_${a}.fasta"
	cd ../
	echo "Round ${a} of polishing complete"
done

echo "Done"

