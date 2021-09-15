#!/bin/bash

# Job name:
#SBATCH --job-name=hisat2
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=18
#
# Memory per node:
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=3:50:00

export PATH=/mnt/home/liujind1/software/hisat2-2.1.0:$PATH
export PATH=/mnt/home/liujind1/software/sraToolkit-2.9.4-1/bin:$PATH
export PATH=/mnt/home/liujind1/software/samtools-1.9/bin:$PATH

export PATH=/mnt/home/liujind1/software/Python-2.7.13/bin:$PATH
export PATH=/mnt/home/liujind1/software/customized:$PATH

export LD_LIBRARY_PATH=/mnt/home/liujind1/software/Python-2.7.13/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/mnt/home/liujind1/software/lapack-3.8.0:$LD_LIBRARY_PATH
export LAPACK=/mnt/home/liujind1/software/lapack-3.8.0
export LD_LIBRARY_PATH=/mnt/home/liujind1/software/BLAS-3.8.0:$LD_LIBRARY_PATH
export export BLAS=/mnt/home/liujind1/software/BLAS-3.8.0/libfblas.a

date > /mnt/scratch/liujind1/fastqDir/log.runTime.SRX1254122
perl /mnt/home/liujind1/software/customized/psi.in.experiment.pl SRX1254122 SRR2432521,SRR2432522 /mnt/scratch/liujind1/fastqDir /mnt/home/liujind1/workAS/01-cattle/database 24
date >> /mnt/scratch/liujind1/fastqDir/log.runTime.SRX1254122

#hisat2 -x /mnt/home/liujind1/workAS/01-cattle/database/genome -1 /mnt/scratch/liujind1/fastqDir/SRR2432521_1.fastq,/mnt/scratch/liujind1/fastqDir/SRR2432522_1.fastq -2 /mnt/scratch/liujind1/fastqDir/SRR2432521_2.fastq,/mnt/scratch/liujind1/fastqDir/SRR2432522_2.fastq --phred33 --rna-strandness RF -p 24 --novel-splicesite-outfile /mnt/scratch/liujind1/fastqDir/SRX1254122.novel.splicesite.rpt.by.hisat.txt  --summary-file /mnt/scratch/liujind1/fastqDir/SRX1254122.summary.align.txt
