#!/bin/bash
#PBS -N Step5b
#PBS -S /bin/bash
#PBS -l walltime=122:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=128gb
#PBS -o /scratch/path/logfiles/log.step5d.out
#PBS -e /scratch/path/logfiles/log.step5d.err
#PBS -d /scratch/path
module load java-jdk/1.8.0_92
module load gatk/3.7
echo START
for input in `ls /scratch/path/*step4.bam`; do java -jar /apps/software/java-jdk-1.8.0_92/gatk/3.7/GenomeAnalysisTK.jar -T PrintReads -R /gpfs/data/godley-lab/WES_analysis/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa -I $input -BQSR ${input%_step4.bam}"_step5a_BaseRecalibrator.table" -o ${input%_step4.bam}"_step5d.bam" ; done
echo END

