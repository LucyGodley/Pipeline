#!/bin/bash
#PBS -N Step5A
#PBS -S /bin/bash
#PBS -l walltime=122:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=128gb
#PBS -o /scratch/path/logfiles/log.step5a.out
#PBS -e /scratch/path/logfiles/log.step5a.err
#PBS -d /scratch/path
module load java-jdk/1.8.0_92
module load gatk/3.7
echo START
for input in `ls /scratch/path/*_step4.bam`; do java -jar /apps/software/java-jdk-1.8.0_92/gatk/3.7/GenomeAnalysisTK.jar -T BaseRecalibrator -R /gpfs/data/godley-lab/WES_analysis/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa -I $input -knownSites /gpfs/data/godley-lab/WES_analysis/dbsnp/All_20170710_sort.vcf -o ${input%_step4.bam}"_step5a_BaseRecalibrator.table" -U ALLOW_SEQ_DICT_INCOMPATIBILITY; done
echo END
