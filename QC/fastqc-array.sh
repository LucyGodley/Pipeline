#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-7
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH -J 21-2490
#SBATCH -o fastqc_%A_%a.out
#SBATCH -e fastqc_%A_%a.err

###################################

#AUTHOR-ASHWIN KOPPAYI
#VERSION-v1
#DATE-03-10-2023

###################################
module load fastqc/0.12.0
mkdir fastqc

ls *.fastq.gz > list_of_files.txt
IFS=$'\n' read -d '' -r -a lines < list_of_files.txt

fastqc -o ./fastqc -t 6 ${lines[$SLURM_ARRAY_TASK_ID]}
