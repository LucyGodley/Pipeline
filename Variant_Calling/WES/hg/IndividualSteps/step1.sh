#!/bin/bash
#PBS -N Step1
#PBS -S /bin/bash
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=256gb
#PBS -o /scratch/path/logfiles/log.step1.out
#PBS -e /scratch/path/logfiles/log.step1.err 
#PBS -d /scratch/path/
module load java-jdk/1.8.0_92
module load picard/2.8.1

for input in `ls /*_R1_001.fastq.gz`; do sample_n=("`echo $input | awk -F"_" '{print $4}'`"); num=("`echo $input | grep  -oP '(?<=L00)\d+(?=\_)'`"); java -Djava.io.tmpdir=${input%_R1_001.fastq.gz}"_tmp" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar FastqToSam FASTQ=$input FASTQ2=${input%_R1_001.fastq.gz}"_R2_001.fastq.gz" OUTPUT=${input%_R1_001.fastq.gz}"ubam.bam" READ_GROUP_NAME=Godley SAMPLE_NAME=$sample_n LIBRARY_NAME=L1 PLATFORM_UNIT=Godley.$num PLATFORM=illumina ;done

echo END
