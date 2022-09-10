#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name download_microbe_database
#SBATCH --output=%x-%j.SLURMout


mkdir microbe-database
cd microbe-database

#download files
echo "Downloading sbt.gtdb-rs202.genomic.k31"
wget -O sbt.gtdb-rs202.genomic.k31.zip https://osf.io/dmsz8/download

#unzip files
echo "Unzipping sbt.gtdb-rs202.genomic.k31.zip"
unzip sbt.gtdb-rs202.genomic.k31.zip

rm unzip sbt.gtdb-rs202.genomic.k31.zip

echo "Down"

