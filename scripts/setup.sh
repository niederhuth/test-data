#!/bin/bash -login
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=6gb
#PBS -N setup

cd $PBS_O_WORKDIR
i=$(pwd | sed s/^.*\\///)
db=$(awk -v FS=',' -v x=$i '{if ($1 == x) print $3}' ../../misc/genomes.csv)
module load bowtie2/2.3.1
module load SAMTools/1.5

mkdir ref
mkdir ref/sequences/ ref/methylCseq ref/annotations
cd ref/sequences

#Download genome
#echo "Downloading $i genome"
#python ../../../../scripts/download_genomes.py $i $db
#rm cookies
#mv $i.gff.gz ../annotations

#Prep genome
echo "Setting up $i"
gunzip $i.fa.gz
samtools faidx $i.fa
cut -f1,2 $i.fa.fai > $i.genome
gunzip ../../../../misc/ChrL.fa.gz -c > tmp
if [ -f ChrC.fa.gz ]
then
	gunzip ChrC.fa.gz
	if [ -f ChrM.fa.gz ]
	then
		gunzip ChrM.fa.gz
		cat $i.fa ChrC.fa ChrM.fa tmp > ../methylCseq/tmp
	else
		cat $i.fa ChrC.fa tmp > ../methylCseq/tmp
	fi
else
	cat $i.fa tmp > ../methylCseq/tmp
fi
rm tmp

#methylCseq index
cd ../methylCseq
python ../../../../scripts/fix_fasta.py -i tmp -o $i.fa
rm tmp
samtools faidx $i.fa
cut -f1,2 $i.fa.fai > $i.genome
methylpy build-reference --input-files $i.fa \
	--output-prefix $i --bowtie2 True

#gzip sequences
cd ../sequences
for x in *fa
do
	gzip $x
done
