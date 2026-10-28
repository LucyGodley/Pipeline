#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=1G

###October 22, 2026
###Author: Paige Kulzer
###
###This script calls nextflow to run a custom whole-exome sequencing pipeline.
###To view a HELP message displaying the list of flags you can use to customize 
###this workflow, add the "--help" flag below.
###e.g. nextflow run main.nf --help

module load multiqc/1.25.1
module load python/3.10.1
module load fastqc/0.12.0
module load nextflow/25.10.3

nextflow run main.nf \
  -profile hpc_slurm \
  --run_dir=/path/to/run \
  --sample_dir=/path/to/run/samples
