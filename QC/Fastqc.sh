#!/bin/bash
module load java-jdk/1.10.0_1
module load fastqc/0.11.7

for files in *fastq.gz;
do 
echo $files
fastqc -t 8 $files
done
