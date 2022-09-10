#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name wtdbg2
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
polish="TRUE"
edge_min=2 #min read coverage for edges, use 2 for low coverage, 4 for high, default is 3

#This should match the dataype in the misc/samples.csv file
#Options include:
#"ont" = Raw Nanopore
#"ont-cor" = Corrected Nanopore
#Need to add PacBio options
datatype="ont" 

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
output="dbg"
path2="wtdbg2"
reads="fastq/${datatype}/clean.fastq.gz"

#Declare reads
echo "reads: ${reads}"

#Set wtdbg2 options based on datatype
if [ ${datatype} = "ont" ]
then
	dt="ont"
	echo "Raw Nanopore reads, using wtdbg2 argument ${dt}"
elif [ ${datatype} = "ont-cor" ]
then
	dt="ont"
	echo "Corrected Nanopore reads, using wtdbg2 argument ${dt}"
else
	echo "Do not recognize ${datatype}"
	echo "Please check and resubmit"
fi

#Get genome size estimate
genomeSize=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $9}' \
	${path1}/samples.csv)

#Run wtdbg2
if [ -d ${path2} ] 
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi
if [ -s ${output}.ctg.lay.gz ]
then
	echo "wtdbg2 output found, proceeding to consensus step"
	echo "To rerun this step, delete ${path2}/${output} and resubmit"
else
	echo "Beginning wtdbg2 assembly"
	wtdbg2 \
		-i ../${reads} \
		-fo ${output} \
		-t ${threads} \
		-x ${dt} \
		-g ${genomeSize} \
		--edge-min ${edge_min} \
		--rescue-low-cov-edges
fi
if [ -s raw/${output}.raw.fa ]
then
	echo "Raw assembly fasta found, proceeding to polishing step"
	echo "To rerun this step, delete ${path2}/raw/${output}.raw.fa and resubmit"
else
	mkdir raw
	echo "Generating raw assembly with wtpoa-cns"
	wtpoa-cns \
	-i ${output}.ctg.lay.gz \
	-fo raw/${output}.raw.fa \
	-t ${threads}
fi
if [ ${polish} == "TRUE" ]
then
	echo "Polishing reads"
	if [ -s polished/${output}.cns.fa ]
	then
		echo "Polished fasta found"
		echo "To rerun this step, delete ${path2}/polished/${output}.cns.fa and resubmit"
	else
		mkdir polished
		if [ -s polished/${output}.bam ]
		then
			echo "Mapped reads found, proceeding to consensus step"
			echo "To rerun this step, delete ${path2}/polished/${output}.bam and resubmit"
		else
			echo "Mapping reads to raw assembly"
			minimap2 \
				-t ${threads} \
				-ax map-${dt} \
				-r2k \
				raw/${output}.raw.fa \
				../${reads} | samtools sort -@4 > polished/${output}.bam
		fi
		echo "Calling consensus with wtpoa-cns"	
		samtools view -F0x900 polished/${output}.bam | \
		wtpoa-cns \
			-t ${threads} \
			-d raw/${output}.raw.fa \
			-i - \
			-fo polished/${output}.cns.fa
		cp polished/${output}.cns.fa ${output}.cns.fa
	fi
else
	cp raw/${output}.raw.fa ./
fi

echo "Done"
