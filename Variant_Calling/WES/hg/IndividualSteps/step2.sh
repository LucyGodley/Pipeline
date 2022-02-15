#!/bin/bash
#PBS -N Step
#PBS -S /bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=200gb
#PBS -o /scratch/path/logfiles/log.step2.out
#PBS -e /scratch/path/logfiles/log.step2.err
#PBS -d /scratch/path
module load java-jdk/1.8.0_92
module load picard/2.8.1
echo START
for input in `ls /scratch/path/*ubam.bam`; do java -Djava.io.tmpdir=${input%ubam.bam}"_tmp_step2" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MarkIlluminaAdapters I=$input O=${input%ubam.bam}"_step2.bam" M=${input%ubam.bam}"_step2_metrics.txt" ;done
echo END
