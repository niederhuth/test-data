#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name flye
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
iterations=2

#This should match the dataype in the misc/samples.csv file
#Options include:
#"ont" = Raw Nanopore
#"ont-cor" = Corrected Nanopore
#"pac" = raw PacBio
#"pac-cor" = Corrected PacBio
#"hifi" = PacBio HiFi
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
path2="flye"
reads="fastq/${datatype}/clean.fastq.gz"

#Declare reads
echo "reads: ${reads}"

#Set flye options based on datatype
if [ ${datatype} = "ont" ]
then
	dt="--nano-raw"
	echo "Raw Nanopore reads, using flye argument ${dt}"
elif [ ${datatype} = "ont-cor" ]
then
	dt="--nano-corr"
	echo "Corrected Nanopore reads, using flye argument ${dt}"
elif [ ${datatype} = "pac" ]
then
	dt="--pacbio-raw"
	echo "Raw PacBio reads, using flye argument ${dt}"
elif [ ${datatype} = "pac-cor" ]
then
	dt="--pacbio-corr"
	echo "Corrected PacBio reads, using flye argument ${dt}"
elif [ ${datatype} = "hifi" ]
then
	dt="--pacbio-hifi"
	echo "PacBio Hifi reads, using flye argument ${dt}"
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

#Run flye
if  [ -s ${path2}/flye.log ]  
then
	echo "Previous flye assembly log detected, restarting from previous run"
	echo "To start from the beginning, please delete the directory ${path2} and resubmit"
	flye \
		--out-dir ${path2} \
		${dt} ${reads} \
		--genome-size ${genomeSize} \
		--iterations ${iterations} \
		--threads ${threads} \
		--resume
else
	echo "Beginning flye assembly"
	flye \
		--out-dir ${path2} \
		${dt} ${reads} \
		--genome-size ${genomeSize} \
		--iterations ${iterations} \
		--threads ${threads} 
fi

echo "Done"
