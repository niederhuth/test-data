#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH --job-name racon
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
rounds=3
datatype="ont"
input="" #Can set to empty and script will find fasta in directory submitted

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)

#Extract reads from assembly job report
reads="$(grep reads: ${path2}/job_reports/${assembly}-*.SLURMout | head -1 | cut -d ' ' -f2)"

#Change preset based on datatype
if [ ${datatype} = "ont" ]
then
	preset="map-ont"
elif [ ${datatype} = "ont-cor" ]
then
	preset="map-ont"
elif [ ${datatype} = "pac" ]
then
	preset="map-pb"
elif [ ${datatype} = "pac-cor" ]
then
	preset="map-pb"
elif [ ${datatype} = "hifi" ]
then
	preset="map-pb"
else
	echo "Do not recognize ${datatype}"
	echo "Please check and resubmit"
fi

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
	path3="racon_${a}"
	echo "Round ${a} of polishing" 
	if [ -d ${path3} ]
	then
		echo "Directory ${path3} exists."
		echo "Checking files."
		cd ${path3}
	else
		mkdir ${path3}
		cd ${path3}
	fi
	#Align data with minimap2
	if [ -s round_${a}.paf ]
	then
		echo "Round ${a} alignment found"
		echo "To rerun this step, delete ${path3}/round_${a}.paf and resubmit"
	else
		echo "Aligning with minimap2"
		minimap2 \
			-t ${threads} \
			-x ${preset} \
			${ref} \
			../../${reads} > round_${a}.paf
	fi
	#Polish with Racon
	if [ -s racon_${a}.fa ]
	then
		echo "Round ${a} polishing found"
		echo "To rerun this step, delete ${path3}/racon_${a}.fa and resubmit"
	else
		echo "Polishing data with Racon"
		racon \
			--include-unpolished \
			-m 8 \
			-x -6 \
			-g -8 \
			-w 500 \
			-t ${threads} \
			../../${reads} \
			round_${a}.paf \
			${ref} > racon_${a}.fa
	fi
	ref="../${path3}/racon_${a}.fa"
	cd ../
	echo "Round ${a} of polishing complete"
done

echo "Done"

