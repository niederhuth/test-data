#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name allmaps
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
distance=rank #cM or rank
mask_regions="$(pwd | sed s/data.*/misc/)/genetic_map/$(pwd | sed s/.*$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)\\/// | sed s/\\/.*//)_mask_regions.bed" #bed file of regions to exclude which may introduce errors, e.g. known inversions
weights= #path to set of alternative weights.txt file
primers=TRUE #Paired primer sequences for genetic markers
primer_sets="Lowry_et_al" #List of primers & associated genetic map
primer_max_dist=5000 #max distance for primers to be separated
synteny=TRUE #Use synteny, right now this assumes same species
gff= #gff file of annotations for synteny, if left blank will look in current directory
synt_ref="IM62 NONTOL TOL" #List of genomes to use for synteny, these are assumed to be the same species
markers=FALSE #Have not implemented this option
optical=FALSE #Have not implemented this option

#quick_synt
#If you do not yet have annotations for your genome, you can use this to map transcript sequences
#From another genome and use this as a set of markers
#I don't really recommend this
quick_synt=FALSE #Use quick_synt TRUE or FALSE
quick_synt_ref="IM62 NONTOL TOL" #List of genomes to use for quick_synt, these are assumed to be of same species
chr_list= #Names of sequences (i.e. chromosome names) to use with quick_synt

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
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3="allmaps"
path4=$(pwd | sed s/${species}\\/.*/${species}/)

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

#Look for gff file if synteny set to TRUE
if [ ${synteny} = "TRUE" ]
then
	if [ -z ${gff} ]
	then
		echo "No input gffprovided, looking for gff"
		if ls *.gff >/dev/null 2>&1
		then
			gff=$(ls *gff | sed s/.*\ //)
			echo "GFF file ${gff} found"
		else
			echo "No gff file found, please check and restart"
		fi
	else
		echo "Input fasta: ${input}"
	fi
fi

#Make and cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Copy input fasta
if [ ! -s input.fa ]
then
	cp ../${input} input.fa
fi

#Aligning primer data with bowtie
if [ ${primers} = "TRUE" ]
then
	#Make bowtie index
	if [ -s bowtie_index/input.rev.2.ebwt ]
	then
		echo "bowtie2 index found"
	else
		echo "Building bowtie2 index"
		mkdir bowtie2_index
		bowtie2-build \
			--quiet \
			--threads ${threads} \
			input.fa \
			bowtie2_index/input
	fi
	for i in ${primer_sets}
	do
		if [ -s ${i}_primers/primers.bed ]
		then
			echo "${i}_primers/primers.bed already exists, skipping"
			echo "To repeat this step, delete ${i}_primers/primers.bed and resubmit"
		else
			mkdir ${i}_primers
			#Align primers
			echo "Aligning primer data with bowtie"
			bowtie2 \
				-p ${threads} \
				--sam-nohead \
				--no-unal \
				--very-sensitive \
				-f \
				-X ${primer_max_dist} \
				-x bowtie2_index/input \
				-1 ${path1}/genetic_map/${i}_forward_primers.fa \
				-2 ${path1}/genetic_map/${i}_reverse_primers.fa \
				-S ${i}_primers/primers.sam

			#Format primers.bed
			for a in $(sed '1d' ${path1}/genetic_map/${i}_genetic_map.csv)
			do
				M=$(echo ${a} | cut -d ',' -f1)
				LG=$(echo ${a} | cut -d ',' -f2)
				GP=$(echo ${a} | cut -d ',' -f3)
				F=$(grep ${M}_F ${i}_primers/primers.sam | cut -f1,3,4)
				Fchr=$(echo ${F} | cut -d ' ' -f2)
				Fpos=$(echo ${F} | cut -d ' ' -f3)
				R=$(grep ${M}_R ${i}_primers/primers.sam | cut -f1,3,4)
				Rchr=$(echo ${R} | cut -d ' ' -f2)
				Rpos=$(echo ${R} | cut -d ' ' -f3)
				if [[ -n ${F} && -n ${R} ]]
				then
					if [ ${Fchr} == ${Rchr} ]
					then
						#echo "Primers ${M} properly paired"
						if [ ${Fpos} -gt ${Rpos} ]
						then
							echo "${Fchr},${Rpos},${Fpos},LG${LG}:${GP}" | \
							tr ',' '\t' >> ${i}_primers/primers.bed
						elif [ ${Fpos} -lt ${Rpos} ]
						then
							echo "${Fchr},${Fpos},${Rpos},LG${LG}:${GP}" | \
							tr ',' '\t' >> ${i}_primers/primers.bed
						fi
					else
						echo "Primers ${M} improperly paired, skipping"
					fi
				elif [[ -n ${F} ]]
				then
					#echo "${M} forward primer only"
					echo "${Fchr},${Fpos},$(expr ${Fpos} + 1),LG${LG}:${GP}" | \
					tr ',' '\t' >> ${i}_primers/primers.bed
				elif [[ -n ${R} ]]
				then
					#echo "${M} reverse primer only"
					echo "${Rchr},${Fpos},$(expr ${Rpos} + 1),LG${LG}:${GP}" | \
					tr ',' '\t' >> ${i}_primers/primers.bed
				else
					echo "Primers ${M} missing data" >> ${i}_primers/missing_primers.txt
				fi
				F=""
				R=""
			done
		fi
		#Filter regions if mask_regions provided
		if [ ! -z ${mask_regions} ]
		then
			echo "Filtering bed file for masked regions"
			bedtools intersect -v \
				-a ${i}_primers/primers.bed \
				-b ${mask_regions} > ${i}_primers/tmp
			mv ${i}_primers/tmp ${i}_primers/primers.bed
		fi
		#Add to list of data for allmaps
		position_data="${i}_primers/primers.bed ${position_data}"
	done
fi

#Synteny 
if [ ${synteny} = "TRUE" ]
then
	#Iterate of references and prep files
	for ref in ${synt_ref}
	do
		if [ -s synt_${ref}/${ref}_synteny.bed ]
		then
			echo "synt_${ref}/${ref}_synteny.bed already exists, skipping"
			echo "To repeat this step, delete synt_${ref}/${ref}_synteny.bed and resubmit"
		else
			echo "Performing synteny analysis for ${ref}"
			mkdir synt_${ref}
			cd synt_${ref}
			echo "Formatting files"
			#Format query gff file
			python -m jcvi.formats.gff bed \
				--primary_only \
				--type=mRNA \
				--key=Name \
				../../${gff} \
				-o input.bed
			#Make cds file
			gffread \
				-x tmp_cds.fa \
				-g ../../${input} \
				../../${gff}
			#Format query cds fasta file
			python -m jcvi.formats.fasta format \
				tmp_cds.fa \
				input.cds
			rm tmp_cds.fa
			#Format target gff file
			python -m jcvi.formats.gff bed \
				--primary_only \
				--type=mRNA \
				--key=Name \
				${path4}/${ref}/ref/annotations/${ref}-v*.gff \
				-o ref.bed
			#Format target cds fasta file
			python -m jcvi.formats.fasta format \
				${path4}/${ref}/ref/annotations/${ref}-v*-cds-primary.fa \
				ref.cds
			#Run synteny analysis
			echo "Running synteny analysis"
			python -m jcvi.compara.catalog ortholog \
				--no_strip_names \
				ref input
			#Build synteny.bed
			echo "Building synteny bed for allmaps"
			python -m jcvi.assembly.syntenypath bed \
				--switch \
				ref.input.anchors \
				-o ${ref}_synteny.bed
			cd ../
		fi
		#Filter regions if mask_regions provided
		if [ ! -z ${mask_regions} ]
		then
			echo "Filtering bed file for masked regions"
			bedtools intersect -v \
				-a synt_${ref}/${ref}_synteny.bed \
				-b ${mask_regions} > synt_${ref}/tmp
			mv synt_${ref}/tmp synt_${ref}/${ref}_synteny.bed
		fi
		#Add to list of data for allmaps
		position_data="synt_${ref}/${ref}_synteny.bed ${position_data}"
	done
fi

#Quick Synteny data
if [ ${quick_synt} = "TRUE" ]
then
	for ref in ${quick_synt_ref}
	do
		if [ -s quick_synt_${ref}/${ref}_synteny.bed ]
		then
			echo "quick_synt_${ref}/${ref}_synteny.bed already exists, skipping"
			echo "To repeat this step, delete quick_synt_${ref}/${ref}_synteny.bed and resubmit"
		else
			#Copy over data
			mkdir quick_synt_${ref}
			cp ${path4}/${ref}/ref/${ref}-v*.fa quick_synt_${ref}/ref.fa
			cp ${path4}/${ref}/ref/annotations/${ref}-v*-cds-primary.fa quick_synt_${ref}/ref-cds-primary.fa
			cp input.fa quick_synt_${ref}/input.fa
			for i in input ref
			do
				#Map transcripts to fasta
				minimap2 \
					-x splice quick_synt_${ref}/${i}.fa \
					quick_synt_${ref}/ref-cds-primary.fa > quick_synt_${ref}/${i}.paf 

				#Convert to bed file
				awk -v OFS="\t" '{print $6,$8,$9,$1,100,$5}' \
				quick_synt_${ref}/${i}.paf > quick_synt_${ref}/${i}_aln.bed

				#Retain only unique alignments
				cut -f4 quick_synt_${ref}/${i}_aln.bed | sort | uniq -c | sed 's/^ *//' | \
				awk -v FS=" " '{if ($1 == 1) print $2}' > quick_synt_${ref}/tmp
				fgrep -f quick_synt_${ref}/tmp quick_synt_${ref}/${i}_aln.bed > quick_synt_${ref}/${i}.bed
				rm quick_synt_${ref}/tmp
			done

			#Filter out only chromosomes in chr_list, make anchor_list
			if [[ ${chr_list} ]]
			then
				echo ${chr_list} | tr ' ' '\n' > quick_synt_${ref}/chr_list.txt
				fgrep -w -f quick_synt_${ref}/chr_list.txt quick_synt_${ref}/ref.bed | \
				cut -f4 > quick_synt_${ref}/anchor_list
			else
				cut -f4 quick_synt_${ref}/ref.bed > quick_synt_${ref}/anchor_list
			fi

			#Make anchors
			fgrep -f quick_synt_${ref}/anchor_list quick_synt_${ref}/input.bed | cut -f4 | \
			awk -v OFS="\t" '{print $1,$1,100}' > quick_synt_${ref}/ref.input.1x1.anchors

			#Build synteny.bed
			python -m jcvi.assembly.syntenypath bed \
				--switch \
				quick_synt_${ref}/ref.input.1x1.anchors \
				-o quick_synt_${ref}/${ref}_quick_synteny.bed
			
			#Cleanup
			rm quick_synt_${ref}/input.fa quick_synt_${ref}/ref-cds-primary.fa quick_synt_${ref}/ref.fa
		fi
		#Filter regions if mask_regions provided
		if [ ! -z ${mask_regions} ]
		then
			echo "Filtering bed file for masked regions"
			bedtools intersect -v \
				-a quick_synt_${ref}/${ref}_quick_synteny.bed \
				-b ${mask_regions} > quick_synt_${ref}/tmp
			mv quick_synt_${ref}/tmp quick_synt_${ref}/${ref}_quick_synteny.bed
		fi
		#Add to list of data for allmaps
		position_data="quick_synt_${ref}/${ref}_quick_synteny.bed ${position_data}"
	done
fi

#Merge files and create weights file
echo "Merging position data files"
python -m jcvi.assembly.allmaps mergebed \
	${position_data} -o allmaps.bed

#if weights.txt is specified, copy it over
if [ -z ${weights} ]
then
	echo "Using allmaps mergebed specified weights.txt"
else
	echo "Using weights file: ${weights}"
	cp ${weights} ./weights.txt
fi

#Run allmaps path
echo "Running allmaps path"
python -m jcvi.assembly.allmaps path \
	--cpus=${threads} \
	--distance=${distance} \
	--noplot \
	allmaps.bed input.fa

#Cleanup
mkdir fasta_other
mv allmaps.chr.fasta fasta_other
mv allmaps.unplaced.fasta fasta_other
mv input.fa fasta_other

echo "Done"




