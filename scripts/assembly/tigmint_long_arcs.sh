#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name tigmint-long-arcs
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20 #doesn't seem to want to use more than 6
cut=500 #cut length for long reads
span=auto #Number of spanning molecules threshold. Set span=auto to automatically select
dist=auto #Max dist between reads to be considered same molecule. auto to automatically calculate
window=1000 #Window size (bp) for checking spanning molecules
minsize=2000 #Minimum molecule size
trim=0 #Number of bases to trim off contigs following cuts
datatype="ont" #ont or pb
input= #input fasta, if left blank, will look for it in current directory, mutually exclusive with input_dir
input_dir="pilon" #common directory name, e.g. pilon or polca to look for assemblies, eclusive with "input"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)

#Get genome size estimate
genomeSize=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $9}' \
	${path1}/samples.csv)
genomeSize2=$(python -c "print(${genomeSize/g/} * 1000000000)")

#Maker sure only one is specified
if [ ${input} ] && [ ${input_dir} ]
then
	echo "Error: input & input_dir cannot both be set!"
	echo "Please set only one"
	exit 1
fi

#Get list of assemblies to work on
if [ -z input_dir ]
then
	assembly_list=${assembly}
else
	assembly_list=${input_dir}*
fi

#Iterate over assembly list and run tigmint-long
for i in ${assembly_list}
do

	#set directory name
	if [ -z ${input_dir} ]
	then
		path3="tigmint_long_arcs"
		path4=".."
	else
		path3="tigmint_long_arcs_${i}"
		path4="../${i}"
	fi

	#Make and cd to output directory
	if [ -d ${path3} ]
	then
		echo "Previous run for $i already found, skipping"
		echo "To rerun this step, delete directory ${path3} and resubmit"
	else
		mkdir ${path3}
		cd ${path3}

		#Look for fasta file, there can only be one!
		if [ -z ${input} ] && [ ${input_dir} ]
		then
			echo "No input fasta provided, looking for fasta"
			if ls ${path4}/*.fa >/dev/null 2>&1
			then
				input=$(ls ${path4}/*fa | sed s/.*\\///)
				name=${input/.fa/}
				echo "Fasta file ${input} found"
			elif ls ${path4}/*.fasta >/dev/null 2>&1
			then
				input=$(ls ${path4}/*fasta | sed s/.*\\///)
				name=${input/.fasta/}
				echo "Fasta file ${input} found"
			elif ls ${path4}/*.fna >/dev/null 2>&1
			then
				input=$(ls ${path4}/*fna | sed s/.*\\///)
				name=${input/.fna/}
				echo "Fasta file ${input} found"
			else
				echo "No fasta file found, please check and restart"
			fi
		else
			echo "Input fasta: ${input}"
		fi

		#Copy and rename files...because of stupid eccentricities of some code
		cp ${path4}/${input} ${name}.fa
		cp ${path2}/fastq/${datatype}/clean.fastq.gz reads.fq.gz

		#Run tigmint
		echo "Running tigmint-long on ${i}"
		tigmint-make tigmint-long arcs \
			draft=${name} \
			reads=reads \
			longmap=${datatype} \
			cut=${cut} \
			span=${span} \
			dist=${dist} \
			window=${window} \
			minsize=${minsize} \
			trim=${trim} \
			G=${genomeSize2/.0/} \
			t=${threads}

		#Clean some stuff up for downstream analyses
		unlink ${name}.cut${cut}.tigmint.fa
		rm ${name}.fa ${name}.fa.fai reads.fq.gz

		#rename fasta & bed file to something more easily handled
		long_part=cut${cut}.molecule.size${minsize}.trim${trim}.window${window}.span${span}.breaktigs
		name2=${name}.reads.${long_part}
		mv ${name2}.fa ${name}_tigmint.fa
		mv ${name2}.bed ${name}_tigmint.bed

		cd ../
		echo "tigmint-long arcs on ${i} complete"

		if [ ${input_dir} ]
		then
			input=
		fi
	fi
done

echo "Done"
