#!/bin/bash
#PBS -N Step3b
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=180gb
#PBS -o /scratch/path/logfiles/log.step3b.out
#PBS -e /scratch/path/logfiles/log.step3b.err
#PBS -d /scratch/path
module load gcc/6.2.0
module load java-jdk/1.8.0_92
module load picard/2.8.1
module load bwa/0.7.15
echo START
for input in `ls /scratch/path/*_3a_forBWA.fastq`; do /apps/software/gcc-6.2.0/bwa/0.7.15/bwa mem -M -t 7 -p /gpfs/data/godley-lab/WES_analysis/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz $input > ${input%_3a_forBWA.fastq}"_3b_bwa_mem.sam" ; done
echo END
