#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30GB
#SBATCH --constraint="NOAUTO:intel18&v100"
#SBATCH --job-name guppy
#SBATCH --output=%x-%j.SLURMout

#Note: Because of restrictions on guppy_basecaller, you need to have your own guppy setup
#This script will only work on the MSU HPCC.

cd ${PBS_O_WORKDIR}
module purge
module use /mnt/home/johnj/software/modulefiles
module load Guppy/4.2.2 

#Set variables
config="dna_r9.4.1_450bps_hac.cfg"
fast5="reads"
out="guppy"

#Run guppy
mkdir ${out}
guppy_basecaller \
	-r \
	-i ${fast5} \
	-s ${out} \
	-c ${config} \
	--num_callers 1 \
	--cpu_threads_per_caller 6 \
	--trim_strategy dna \
	--calib_detect \
	--qscore_filtering \
	--min_qscore 7 \
	-x auto --gpu_runners_per_device 12 \
	--compress_fastq \
	--disable_pings

echo "Done"