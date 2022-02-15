#!/bin/bash
#PBS -N Step4
#PBS -S /bin/bash
#PBS -l walltime=82:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=80gb
#PBS -o /scratch/path/logfiles/log.step4.out
#PBS -e /scratch/path/logfiles/log.step4.err
#PBS -d /scratch/path
module load java-jdk/1.8.0_92
module load picard/2.8.1
echo START
for input in `ls /scratch/path/*_mapped.bam`; do java -Djava.io.tmpdir=${input%_L004_mapped.bam}"_tmp_step4" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MarkDuplicates INPUT=$input OUTPUT=${input%_L004_mapped.bam}"_step4.bam" METRICS_FILE=${input%_L004_mapped.bam}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true ; done
echo END
