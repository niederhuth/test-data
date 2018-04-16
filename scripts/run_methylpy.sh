#!/bin/bash -login
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=6
#PBS -l mem=100gb
#PBS -N methylpy

cd $PBS_O_WORKDIR

sample=$(pwd | sed s/^.*\\///)
species=$(echo $sample | sed s/-.*//)
refChr=$(awk -v FS=',' -v x=$sample '{if ($1 == x) print $4}' ../../misc/genomes.csv)
echo $refChr

module load SAMTools/1.5
module load bowtie2/2.3.1
module load Java/1.8.0_31
module load SRAToolkit/2.8.2

echo "Running $sample"
mkdir methylCseq
mkdir methylCseq/fastq
cd methylCseq/fastq
if ls *fastq.gz >/dev/null 2>&1
then
	echo "Data present"
else
	echo "Downloading from SRA"
	python ../../../../scripts/download_fastq.py "$sample"
fi
if ls *sra >/dev/null 2>&1
then
	for i in *sra
	do
		fastq-dump --split-3 $i
		rm $i
	done
fi
echo "Unpacking fastqs"
for i in *fq.gz *fastq.gz
do
	gunzip $i
done
for i in *fq
do
	name=$(echo $i | sed s/.fq$/.fastq/)
	mv $i $name
done
cd ../

if ls fastq/*_2.fastq >/dev/null 2>&1
then
	echo "Data is paired-end"
	echo "Running methylpy"
	methylpy paired-end-pipeline \
	--read1-files fastq/*_1.fastq \
	--read2-files fastq/*_2.fastq \
	--libraries "libA" \
	--sample $sample \
	--forward-ref ../../$species/ref/methylCseq/"$species"_f \
	--reverse-ref ../../$species/ref/methylCseq/"$species"_r \
	--ref-fasta ../../$species/ref/methylCseq/$species.fa \
	--num-procs 6 \
	--sort-mem 9 \
	--trim-reads True \
	--path-to-cutadapt "" \
	--adapter-seq-read1 AGATCGGAAGAGCACACGTCTGAAC \
	--adapter-seq-read2 AGATCGGAAGAGCGTCGTGTAGGGA \
	--max-adapter-removal 1 \
	--overlap-length 3 \
	--error-rate 0.1 \
	--min-qual-score 10 \
	--min-read-len 30 \
	--bowtie2 True \
	--path-to-aligner "" \
	--aligner-options "" \
	--merge-by-max-mapq True \
	--path-to-samtools "" \
	--remove-clonal True \
	--keep-clonal-stats True \
	--path-to-picard /mnt/home/niederhu/.local/java/lib \
	--java-options "" \
	--binom-test True \
	--num-upstream-bases 0 \
	--num-downstream-bases 2 \
	--unmethylated-control $refChr \
	--min-base-quality 1 \
	--min-cov 3 \
	--sig-cutoff .01 \
	--path-to-output "" \
	--keep-temp-files False \
	--generate-allc-file True \
	--generate-mpileup-file True \
	--remove-chr-prefix False \
	--compress-output True
else
	echo "Data is single-end"
	echo "Running methylpy"
	methylpy single-end-pipeline \
	--read-files fastq/*.fastq \
	--libraries "libA" \
	--sample $sample \
	--forward-ref ../../$species/ref/methylCseq/"$species"_f \
	--reverse-ref ../../$species/ref/methylCseq/"$species"_r \
	--ref-fasta ../../$species/ref/methylCseq/$species.fa \
	--num-procs 6 \
	--sort-mem 9 \
	--trim-reads True \
	--path-to-cutadapt "" \
	--adapter-seq AGATCGGAAGAGCACACGTCTG \
	--max-adapter-removal 1 \
	--overlap-length 3 \
	--error-rate 0.1 \
	--min-qual-score 10 \
	--min-read-len 30 \
	--bowtie2 True \
	--path-to-aligner "" \
	--aligner-options "" \
	--merge-by-max-mapq True \
	--path-to-samtools "" \
	--remove-clonal True \
	--keep-clonal-stats True \
	--path-to-picard /mnt/home/niederhu/.local/java/lib \
	--java-options "" \
	--binom-test True \
	--num-upstream-bases 0 \
	--num-downstream-bases 2 \
	--unmethylated-control $refChr \
	--min-base-quality 1 \
	--min-cov 3 \
	--sig-cutoff .01 \
	--path-to-output "" \
	--keep-temp-files False \
	--generate-allc-file True \
	--generate-mpileup-file True \
	--remove-chr-prefix False \
	--compress-output True
fi

rm *mpileup_output.tsv *_reads_no_clonal*.bam* *_libA.metric

echo "Compressing fastqs"
cd fastq
for i in *fastq
do
	gzip $i
done
