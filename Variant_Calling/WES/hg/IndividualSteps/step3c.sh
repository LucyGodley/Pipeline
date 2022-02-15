#!/bin/bash
#PBS -N Step3C
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=128gb
#PBS -o /scratch/path/logfiles/log.step3c.out
#PBS -e /scratch/path/logfiles/log.step3c.err
#PBS -d /scratch/path
module load gcc/6.2.0
module load java-jdk/1.8.0_92
module load picard/2.8.1
module load bwa/0.7.15
echo START
for input in `ls /scratch/path/*_3b_bwa_mem.sam`; do java -Djava.io.tmpdir=${input%_3b_bwa_mem.sam}"_tmp_step3c_BWA" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MergeBamAlignment R=/gpfs/data/godley-lab/WES_analysis/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz ALIGNED_BAM=$input UNMAPPED_BAM=${input%_3b_bwa_mem.sam}"ubam.bam" OUTPUT=${input%_3b_bwa_mem.sam}"_mapped.bam" CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS ; done
echo END
