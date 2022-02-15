#!/bin/bash
#PBS -N Step3A
#PBS -S /bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=100gb
#PBS -o /scratch/path/logfiles/log.step3a.out
#PBS -e /scratch/path/logfiles/log.step3a.err
#PBS -d /scratch/path
module load java-jdk/1.8.0_92
module load picard/2.8.1
echo START
for input in `ls /scratch/path/*_step2.bam`; do java -Djava.io.tmpdir=${input%_step2.bam}"_tmp_step3a_BWA" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar SamToFastq I=$input FASTQ=${input%_step2.bam}"_3a_forBWA.fastq" CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true; done
echo END
