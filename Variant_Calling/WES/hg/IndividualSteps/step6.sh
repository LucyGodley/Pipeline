#!/bin/bash
#PBS -N Step6
#PBS -S /bin/bash
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=100gb
#PBS -o /scratch/path/logfiles/log.step6.out
#PBS -e /scratch/path/logfiles/log.step6.err
#PBS -d /scratch/path
module load java-jdk/1.8.0_92
module load gatk/3.7
echo START
for input in `ls /scratch/path/*step5d.bam`; do java -jar /apps/software/java-jdk-1.8.0_92/gatk/3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R /gpfs/data/godley-lab/WES_analysis/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa -I $input --genotyping_mode DISCOVERY -o ${input%_step5d.bam}".vcf"; done
echo END

