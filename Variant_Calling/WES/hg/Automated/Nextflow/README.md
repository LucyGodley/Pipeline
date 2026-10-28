# Nextflow WES/WGS Analysis Pipeline

## Overview

This Nextflow pipeline takes FASTQ files as input from either WES or WGS experiments and automates the alignment and filtering of that data for rare germline variant discovery. It is designed to run on an HPC cluster using the SLURM resource manager and generates Excel output files for two gene panels:

1. 378-Gene Panel
2. ACMG Secondary Findings (ACMG SF)

This pipeline has been configured such that all *CHEK2* and *DDX41* variants at or below an allele frequency threshold of 5% are reported in the Excel output, while all other variants at all other genes are reported when allele frequency is at or below 1%.

The pipeline supports configurable execution settings. It also generates a Nextflow execution report to help monitor pipeline success, failures, and resource usage. For more information about the specifics of this pipeline, please see 'WES_pipeline.docx'.

## Requirements

Before running the pipeline, ensure the following are available:

- Nextflow installed and accessible from your environment.
- Access to an HPC cluster with SLURM, such as Quest.
- Python and any dependencies required by the gene filtering script.
- All required pipeline scripts, configuration files, and gene panel files.
- Sample folders containing the input data expected by the pipeline.
- Access to a scratch or project directory for sample data and temporary workflow files.

## Directory Structure

All pipeline scripts, configuration files, and supporting files must be located in the same directory as `main.nf`. Sample folders should be located in a scratch or project directory.

An example directory structure is shown below with all required files in the home/ directory:

- home/
  - main.nf
  - nextflow.config
  - submit_nextflow.sh
  -  Gene_filter_v4.py
  - ACMG_SF_v3.3_Gene_list.txt
  - NW_378Genes_Panel_MOD.txt
  - additional pipeline files...

- scratch/
  - run_directory/
    - sample_directory/
      - Sample_001/
      - Sample_002/
      - Sample_003/

#### **Important**: Sample folders are expected to begin with "Sample_" or another prefix separated from the rest of the folder name by an underscore. This naming convention can be modified in the workflow section at the bottom of main.nf.

## Usage

The typical command for running the pipeline is:

`nextflow run main.nf \
    -profile hpc_slurm \
    --run_dir=/path/to/topdir \
    --sample_dir=/path/to/topdir/sampledir`

Replace the example paths with the appropriate directories for your run.

To submit the pipeline as a SLURM job on Quest, run:

`sbatch submit_nextflow.sh`

The `submit_nextflow.sh` script should contain the appropriate SLURM job settings and the Nextflow command. Job-specific SLURM settings can be modified in `nextflow.config`.

Review the submission script and configuration file before launching a production run to ensure the paths and resource allocations are correct.

## Parameters

### Mandatory arguments
| Argument           | Description                                                                      |
|--------------------|----------------------------------------------------------------------------------|
| -profile | Enables the SLURM execution profile for running the pipeline on the HPC cluster. |
| --run_dir          | Path to the top-level directory containing the sample directory.                 |
| --sample_dir       | Path to the directory containing the sample folders to be processed.             |

### Optional arguments
| Argument             | Default                            | Description                                                                                          |
|--------------------------|------------------------------------|------------------------------------------------------------------------|
| --date               | Today's date in yyyyMMdd format    | Date associated with the pipeline run.                                    |
| --gene_filter_script | Gene_filter_v4.py     | Python script used for gene filtering.                                                               |
| --gene_panel_file_1  | ACMG_SF_v3.3_Gene_list.txt         | Gene list used for the ACMG SF panel.                                                                |
| --gene_panel_file_2  | NW_378Genes_Panel_MOD.txt | Gene list used for the 378-gene panel.                                                               |
| --keep_intronic      | Yes                                | Controls whether intronic variants are retained, according to the gene filter script. |
| --folders            | 1                                  | Number of folders containing annotated VCF output, according to the gene filter script.                                     |
| --output_filename_1  | ACMGGenes_Panel_Output.xlsx        | Filename for the ACMG SF panel output.                                                               |
| --output_filename_2  | 378Genes_Panel_Output.xlsx         | Filename for the 378-gene panel output.                                                              |

**Note**: The default values above reflect the values specified in the pipeline's help function. Confirm that they match the current parameter definitions in `nextflow.config`.

Example with optional parameters:
`nextflow run main.nf \
    -profile hpc_slurm \
    --run_dir=/scratch/path/to/run \
    --sample_dir=/scratch/path/to/run/samples \
    --keep_intronic=No \
    --output_filename_1=Filtered_Vars_ACMG_Output.xlsx \
    --output_filename_2=Filtered_Vars_378_Output.xlsx`


## Resuming a Previous Run

Nextflow supports resuming previous runs, allowing completed tasks to be reused when their cached results are available.

To list previous runs, execute:

`nextflow log`

Identify the run name or identifier associated with the run you want to resume, then add `-resume` to the original command:

`nextflow run main.nf \
    -profile hpc_slurm \
    --run_dir=/path/to/topdir \
    --sample_dir=/path/to/topdir/sampledir \
    -resume <run-id>`


If necessary, specify a particular run identifier after `-resume`, following the syntax supported by your Nextflow version. By default, Nextflow will resume the previous run at the beginning of the step it exited on.

Resuming works best when the relevant work directory and cached task outputs remain available. More information about resuming jobs can be found in the offical Nextflow documentation.

## Work Directory and Configuration

Nextflow uses a work directory to stage input files, store intermediate files, and manage task execution.

The work directory can be configured through the `workDir` variable in `nextflow.config`. It should point to a scratch location with sufficient storage and appropriate access permissions.

Review `nextflow.config` before running the pipeline to configure:

- The work directory location.
- The output location for the Nextflow execution report.
- CPU and memory allocations for each task in the pipeline.

Using scratch space for intermediate files helps keep temporary workflow data separate from the pipeline scripts and input files. More information about setting a work directory can be found in the offical Nextflow documentation.

## Outputs

The pipeline generates several output files organized into subdirectories within the run directory:

| Output folder      | Contents                                                         |
|----------------|------------------------------------------------------------------|
| FASTQ/         | Concatenated FASTQ files.              |
| Logfiles/      | FastQC log files.                                                |
| multiqc/       | MultiQC report summarizing quality control results.              |
| Annotated_VCF/ | VCF files annotated with Annovar. |
| GVCFs/         | GVCF files.               |
| Step5b_bam/    | BAM files. |
| VCF/           | VCF files. |

The gene filter script built into the pipeline generates two Excel workbooks, one for each gene panel used:

2. ACMGGenes_Panel_Output.xlsx (ACMG SF gene panel)
2. 378Genes_Panel_Output.xlsx (378-gene panel)

Output filenames for these workbooks can be changed using the `--output_filename_1` and `--output_filename_2` flags. See 'Parameters' above.

## Execution Reports and Logs

Nextflow generates execution information that can be used to assess pipeline performance and identify failed tasks.

The execution report provides information about:

- Overall workflow execution.
- Task success and failure.
- Resource usage by individual tasks.

The report output location can be configured in `nextflow.config`. Review the report and relevant SLURM logs after each run, especially if the pipeline fails or some samples are not processed as expected.

## Troubleshooting
#### Pipeline cannot find sample folders
- Confirm that `--sample_dir` points to the correct directory containing sample folders with input FASTQ files.
- Check that sample folder names begin with Sample_ or follow the configured prefix-plus-underscore convention.
- Review the folder-matching logic in the workflow section of `main.nf`.

#### Gene filtering script or panel files cannot be found
- Confirm that the required files are present in the same directory as `main.nf`, as expected by the pipeline.
- Check the filenames passed to `--gene_filter_script`, `--gene_panel_file_1`, and `--gene_panel_file_2`.
- Verify that the scripts and files have appropriate read permissions.

#### SLURM job fails or remains pending
- Review the SLURM job status and job output/error logs.
- Confirm that the hpc_slurm profile is configured correctly.
- Check the job-specific settings in `submit_nextflow.sh` and task-specific settings in `nextflow.config`.

#### Work directory or temporary file issues
- Ensure sufficient scratch storage is available.
- Verify that the job has permission to read and write in the configured directories.
- When resuming a run, ensure the previous work directory and cached task outputs have not been removed.

#### Pipeline fails during execution
- Review the Nextflow execution report.
- Inspect the intermediate files and logs in the specified work directory.
- Check the SLURM logs for resource-limit errors or job termination.

## Best Practices
- Verify input paths and sample folder names before launching a run.
- Keep the pipeline scripts and supporting gene panel files together.
- Use scratch space for sample inputs and intermediate workflow files as required by the pipeline.
- Review the config and SLURM resource settings before processing large batches.
- Retain execution reports and relevant logs for troubleshooting and reproducibility.