#! /usr/bin/env nextflow

if (params.help) {
  log.info """
      Usage:
      The typical command for running the pipeline is as follows:
      nextflow run main.nf -profile hpc_slurm --run_dir=/path/to/topdir --sample_dir=/path/to/topdir/sampledir
      
      To submit a nextflow job on Quest, use:
      sbatch submit_nextflow.sh
      
      Mandatory arguments:
      -profile                Allows nf to use the SLURM resource manager
      --run_dir               Path to the top-level directory where your sample directory lives
      --sample_dir            Path to the directory containing sample folders
      
      Optional arguments:
      --date                  Default: Today's date in 'yyyyMMdd' format
      --gene_filter_script    Default: 'Gene_filter_v4.py'
      --gene_panel_file_1     Default: 'ACMG_SF_v3.3_Gene_list.txt'
      --gene_panel_file_2     Default: 'NW_378Genes_Panel_MOD.txt'
      --keep_intronic         Default: 'Yes'
      --folders               Default: 1
      --output_filename_1     Default: 'ACMGGenes_Panel_Output.xlsx'
      --output_filename_2     Default: '378Genes_Panel_Output.xlsx'
      --workDir               Path to the scratch directory where temp files can go
      -with-report            Path to the directory where a run report will be generated
      -resume {run name}      Add this flag to resume a previous nextflow run. Use "nextflow log" on the terminal to list the last run names.
      
      Nextflow utilizes a directory called the work directory to stage files while the pipeline is running. 
      The work directory can be set via the workDir variable in the config file, "nextflow.config". 
      Job-specific SBATCH settings can also be modified in the config file.
      
      Nextflow will output a report at the end of each run which contains more information about the success
      or failure of a run, as well as the resources used by each task. The output directory for this report
      can be specified in the config file.
      
      All files must be in the same directory as the nextflow script (main.nf) with the 
      exception of sample folders which should be in scratch. PLEASE NOTE that sample folders 
      are expected to begin with the string "Sample_" or some other prefix separated by an underscore. 
      This can be modified in the workflow section at the bottom of this script.
  """
  exit 0
  
}

process concatFastqs {
  input:
      tuple val(sample_id), val(dir_names), path(reads)
  
  output:
      tuple val("${sample_id}"), path("${sample_id}*/*.fastq.gz"), emit: combined_fqs

  publishDir(
      path: "${params.run_dir}/output/FASTQ", // Target directory
      mode: 'copy' // Copy files to the target dir
  )
    
  script:
  """
  declare -A FILES           # key: meta_key|R1 or meta_key|R2 -> list of fq paths
  declare -A SNUMS           # key: meta_key|Rx -> list of S numbers (numeric)
  declare -A LANES           # key: meta_key|Rx -> list of lane numbers (numeric)
  
  if [[ ${dir_names.size()} -eq 1 ]]; then
        outdir="${sample_id}"
        echo "Sample: ${sample_id}"
        echo "Only one directory."
        echo "Copying FASTQs without modification."

        mkdir -p "\$outdir"
        cp ${reads} "\$outdir/"
        
    else
        outdir="${sample_id}_combined"
        mkdir -p "\$outdir"

        echo "Sample: ${sample_id}"
        echo "Directories: ${dir_names}"
        echo "Number of directories: ${dir_names.size()}"
        
        # 1) Scan all fastqs and fill maps  
        for fq in ${reads}; do
          base=\$fq       # e.g. CML-007_1_223_354_S63_L002_R1_001.fastq.gz
          
          # Split filename on "_"
          IFS='_' read -r -a parts <<< "\$base"
          n=\${#parts[@]}
      
          s_field=\${parts[n-4]}      # e.g. S63
          l_field=\${parts[n-3]}      # e.g. L002
          r_field=\${parts[n-2]}      # e.g. R1
      
          # numeric values
          s_num=\${s_field#S}
          l_num=\${l_field#L}
      
          # group files by (meta_key, lane, read)
          key="${sample_id}|\${l_field}|\${r_field}"
      
          FILES[\$key]+="\$fq "
          SNUMS[\$key]+="\$s_num "
          LANES[\$key]+="\$l_num "
      
        done
      
      # 2) Concat samples
        for key in "\${!FILES[@]}"; do
          echo "key is \$key"
          
          # key here is ID|Lx|Rx; split
          IFS='|' read -r sample lane read_field <<< "\$key"
      
          # gather S numbers
          min_s=""
          for s in \${SNUMS[\$key]}; do
            if [[ -z "\$min_s" || \$s -lt \$min_s ]]; then
              min_s="\$s"
            fi
          done
          
          printf -v s_out "S%d" "\$min_s"
      
          outfile="\${outdir}/${sample_id}_\${s_out}_\${lane}_\${read_field}_001.fastq.gz"
      
          echo "Creating: "
          echo "\$outfile"
      
          printf "%s\n" \${FILES[\$key]} | sort -V > files.tmp
      
          echo "Input files:"
          cat files.tmp
      
          xargs cat < files.tmp > "\$outfile"
        done
      fi
  echo "Done. Concatenated R1/R2 FASTQs only for samples with multiple directories."
  
  """
}

process fastqcArray {
  input:
      tuple val(sample_id), path(reads)
     
  output:
      path "*_fastqc.zip", emit: fastqc_zips
      path "*.err"
      path "*.out"
      
  publishDir (path: "${params.run_dir}/output/Logfiles", mode: 'move', pattern: "*.{err,out}")      
  
  script:
  """
  fastqc -t 6 ${reads}
  mv *.err fastqc_${sample_id}.err
  mv *.out fastqc_${sample_id}.out
  """
}

process multiQC {
  input:
      path "*"
      
  output:
      path "multiqc_*"
      
  publishDir (path: "${params.run_dir}/output/multiqc", mode: 'move', pattern: "*.html")

  script:
  """
  multiqc . --filename multiqc_${params.date}
  """
}

process WESanalysis {
  input:
      tuple val(sample_id), path(reads)
      
  output:
      path "*multianno.txt", emit: multianno
      path "*multianno.vcf"
      path "*g.vcf.*"
      path "*step5b*"
      path "*mapped.bam.vcf*"
      
  publishDir (path: "${params.run_dir}/output/Annotated_VCF", mode: 'copy', pattern: "*multianno*")
  publishDir (path: "${params.run_dir}/output/GVCFs", mode: 'copy', pattern: "*.g.vcf.*")
  publishDir (path: "${params.run_dir}/output/Step5b_bam", mode: 'copy', pattern: "*step5b*")
  publishDir (path: "${params.run_dir}/output/VCF", mode: 'copy', pattern: "*mapped.bam.vcf*")

  script:
  """
  
  ref_path="/projects/p32126/WES_analysis"
  reference="/reference_hg_38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  dbsnp="/dbsnp/hg38/ALL_20141222.dbSNP142_human_GRCh38_sort.snps.vcf"
  
  echo START
  echo "WES Analysis Started" ;
  
  lane=\$(ls *R1*|wc -l);
  
  echo "sample_id is ${sample_id}"
  echo "lane is \$lane"   
  
  echo "Running the WES pipeline" ;
  

                                                                            #STEP1
  for r1 in *_R1_001.fastq.gz; do
    r2=\${r1%_R1_001.fastq.gz}"_R2_001.fastq.gz"

    echo "Processing R1: \$r1"
    echo "Processing R2: \$r2"     
                                                               
    module load picard/2.6.0
    PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"
    
    output1=\${r1%_R1_001.fastq.gz}"_ubam.bam"
    echo "Checking output file=\$output1 for step 1" ;
    if test -f "\$output1"; then
      echo "\$output1 Already exists" ;
      echo "Skipping Step 1 for ${sample_id}" ;
    else
    echo "Step 1 of Analysis Started for ${sample_id}" 
    platform_unit="\$(zless "\$r1" | head -n1 | awk -F ':' '{print \$3 "." \$4}')";
    java -Djava.io.tmpdir=\${r1%_R1_001.fastq.gz}"_001_tmp" -jar \${PICARD} FastqToSam FASTQ=\$r1 FASTQ2=\$r2 OUTPUT=\$output1 READ_GROUP_NAME=UKB SAMPLE_NAME=${sample_id} LIBRARY_NAME=\$(echo \$r1 | awk -F '_' '{print \$(NF-2)}') PLATFORM_UNIT=\$platform_unit PLATFORM=ILLUMINA;
    
      if test -f "\$output1"; then
        echo "\$output1 exists" ;
        echo "Step 1 of Analysis Completed for ${sample_id}" ;
      else
         echo "Error! \$output1 not found" ;
         continue; 
      fi
    fi
  
  
                                                                         #STEP2
  
    module load picard/2.6.0
    PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"
    
    
    output2=\${r1%_R1_001.fastq.gz}"_step2.bam"
    echo "Checking for output file=\$output2 for step 2" ;
    
    if test -f "\$output2"; then
    echo "\$output2 Already exists" ;
    echo "Skipping Step 2 for ${sample_id}" ;                                    else
    echo "Step 2 of Analysis Started for ${sample_id}" ;
    java -Djava.io.tmpdir=\${r1%_R1_001.fastq.gz}"_tmp_step2" -jar \${PICARD} MarkIlluminaAdapters I=\$output1 O=\$output2 M=\${r1%_R1_001.fastq.gz}"_step2_metrics.txt" ;
    output2_metrics=\${output1%_ubam.bam}"_step2_metrics.txt"
      if test -f "\$output2"; then
         echo "\$output2 exists" ;
         echo "Step 2 of Analysis Completed for ${sample_id}" ;
      else
        echo "Error! \$output2 not found!" ;
        continue;
      fi
    fi
   
                                                                       #STEP3A
  
  
    module load picard/2.6.0
    PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"
    
    output3a=\${r1%_R1_001.fastq.gz}"_3a_forBWA.fastq"
    echo "Checking for output file=\$output3a for step 3a" ; 
    
    if test -f "\$output3a"; then 
    echo "\$output3a Already exists" ;
    echo "Skipping Step3a for ${sample_id}" ;
    else 
    echo "Step 3 of Analysis Started for ${sample_id}" ;
    java -Djava.io.tmpdir=\${r1%_R1_001.fastq.gz}"_tmp_step3a_BWA" -jar \${PICARD} SamToFastq I=\$output2 FASTQ=\$output3a CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true;
    
      if test -f "\$output3a"; then
         echo "\$output3a exists" ;
         echo "Step 3a of Analysis Completed for ${sample_id}" ;
      else
         echo"Error! \$output3a not found!";
         continue;
      fi
    fi
  
                                                                     #STEP3B
  
    module load bwa/0.7.17
    BWA="/software/bwa/gcc/0.7.17/bwa"
    output3b=\${r1%_R1_001.fastq.gz}"_3b_bwa_mem.sam";
    echo "Checking for output file=\$output3b for step 3b" ; 
    
    if test -f "\$output3b" ;then
    echo "\$output3b Already exists" ;
    echo "Skipping Step3b for ${sample_id}" ;
    else
    echo "Step 3b of Analysis Started for ${sample_id}" ;
    \${BWA} mem -M -t 7 -p \${ref_path}"\$reference" \$output3a > \$output3b ;
      if test -f "\$output3b"; then
        echo "\$output3b exists" ;
        echo "Step 3b of Analysis Completed for ${sample_id}" ;
      else
        echo "Error!\$output3b not found!";
        continue;
      fi
    fi
  
                                                                   #STEP3C    
  
    module load picard/2.6.0
    PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"
    
    
    
    output3c=\${r1%_R1_001.fastq.gz}"_mapped.bam"
    echo "Checking for output file=\$output3c for step 3c" ;
    
    if test -f "\$output3c" ;then
    echo "\$output3c Already exists" ;
    echo "Skipping Step3c for ${sample_id}" ;
    else
    echo "Step 3c of Analysis Started for ${sample_id}" ; 
    java -Djava.io.tmpdir=\${r1%_R1_001.fastq.gz}"_tmp_step3c_BWA" -jar \${PICARD} MergeBamAlignment R=\${ref_path}"\$reference" ALIGNED_BAM=\$output3b UNMAPPED_BAM=\$output1 OUTPUT=\$output3c CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS;
    
    
      if test -f "\$output3c"; then
        echo "\$output3c exists" ;
        echo "Step3c of Analysis Completed for ${sample_id}" ;
      else
        echo "Error! \$output3c not found!" ;
        continue;
      fi
    fi
  done
                                                                  #STEP4 
  module load picard/2.6.0
  PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"
  
  output4=${sample_id}*_step4.bam;
  echo "Checking for output file=\$output4 for step 4"  ;
  if test -f "\$output4" ;then
  echo "\$output4 Already exists" ;
  echo "Skipping Step4 for ${sample_id}" ;
  else
  echo "Step 4 of Analysis Started for ${sample_id}" ;
  case \$lane in
  1) echo "${sample_id} has 1 Lanes" ;
  files=(*mapped.bam)
  echo Number of Mapped Bam files= "\${#files[@]}";
  java -Djava.io/.tmpdir=\${output3c%_mapped.bam}"_tmp_step4" -jar \${PICARD} MarkDuplicates INPUT=\$output3c OUTPUT="\${output3c%_mapped.bam}_step4.bam" METRICS_FILE="\${output3c%_mapped.bam}_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true; 
  output4="\${output3c%_mapped.bam}_step4.bam";   
  ;; 
  2) echo "${sample_id} has 2 Lanes"  ;
  files=(*mapped.bam)
  if [ "\${#files[@]}" -eq "\$lane" ]; then
  echo Number of Mapped Bam files= "\${#files[@]}";
  java -Djava.io/.tmpdir=\${output3c%_mapped.bam}"_tmp_step4" -jar \${PICARD} MarkDuplicates INPUT=\${files[0]} INPUT=\${files[1]} OUTPUT="\${output3c%_mapped.bam}_step4.bam" METRICS_FILE="\${output3c%_mapped.bam}_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true;
  output4="\${output3c%_mapped.bam}_step4.bam"  
  else
  continue 
  fi
  ;;
  3) echo "${sample_id} has 3 Lanes" 
  files=(*mapped.bam)
  if [ "\${#files[@]}" -eq "\$lane" ]; then
  echo Number of Mapped Bam files= "\${#files[@]}";
  java -Djava.io/.tmpdir=\${output3c%_mapped.bam}"_tmp_step4" -jar \${PICARD} MarkDuplicates INPUT=\${files[0]} INPUT=\${files[1]} INPUT=\${files[2]} OUTPUT="\${output3c%_mapped.bam}_step4.bam" METRICS_FILE="\${output3c%_mapped.bam}_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true;
  output4="\${output3c%_mapped.bam}_step4.bam";
  else
  continue;
  fi   
  ;;
  4) echo "${sample_id} has 4 Lanes" 
  files=(*mapped.bam)
  if [ "\${#files[@]}" -eq "\$lane" ]; then
  echo Number of Mapped Bam files= "\${#files[@]}";
  java -Djava.io/.tmpdir=\${output3c%_mapped.bam}"_tmp_step4" -jar \${PICARD} MarkDuplicates INPUT=\${files[0]} INPUT=\${files[1]} INPUT=\${files[2]} INPUT=\${files[3]} OUTPUT="\${output3c%_mapped.bam}_step4.bam" METRICS_FILE="\${output3c%_mapped.bam}_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true; 
  output4="\${output3c%_mapped.bam}_step4.bam";
  else
  continue
  fi
  ;;
  *)
  echo "Flop"
  ;;
  esac
  
    if test -f "\$output4"; then
      echo "\$output4 exists" ;
      echo "Step4 of Analysis Completed for ${sample_id}" ;
    else
      echo "Error! \$output4 not found";
      continue;
    fi
  fi
  
                                                              #STEP5a
  
  module load java/jdk-17.0.2+8
  module load gatk/4.4.0.0
  GATK="/software/gatk/4.4.0.0/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
  
  output5a=\${output4%_step4.bam}"_step5a_BaseRecalibrator.table";
  echo "Checking for output file=\$output5a for step 5a" ;
  
  if test -f "\$output5a" ;then
  echo "\$output5a Already exists" ;
  echo "Skipping Step5a for ${sample_id}" ;
  else
  echo "Step 5a of Analysis for ${sample_id} Started" ;
  java -jar \${GATK} BaseRecalibrator -R \${ref_path}"\$reference" -I \$output4 --known-sites \${ref_path}"\$dbsnp" -O \$output5a
    if test -f "\$output5a"; then
      echo "\$output5a exists" ;
      echo "Step 5a of Analysis Completed" ;
    else
      echo "Error!\$output5a not found!";
      continue
    fi
  fi
  
  
                                                           #STEP5b
  
  module load java/jdk-17.0.2+8
  module load gatk/4.4.0.0
  GATK="/software/gatk/4.4.0.0/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
  
  output5b=\${output4%_step4.bam}"_step5b.bam";
  echo "Checking for output file=\$output5b of step 5b" ; 
  
  if test -f "\$output5b" ;then
  echo "\$output5b Already exists" ;
  echo "Skipping Step5b for ${sample_id}" ;                                                                           else 
  echo "Step 5b of Analysis Started" ;
  java -jar \${GATK} ApplyBQSR -R \${ref_path}"\$reference" -I \$output4 -bqsr \$output5a -O \${output4%_step4.bam}"_step5b.bam" ;
  
    if test -f "\$output5b"; then
      echo "\$output5b exists" ;
      echo "Step 5b of Analysis Completed for ${sample_id}" ;
    else
     echo "Error! \$output5b not found! ";
     continue;
    fi
  fi
   
  
                                                          #STEP6
  
  module load java/jdk-17.0.2+8
  module load gatk/4.4.0.0
  GATK="/software/gatk/4.4.0.0/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
  
  output6=\${output5b%_step5b.bam}".g.vcf.gz"
  echo "Checking for output file=\$output6 for step 6" ; 
  
  if test -f "\$output6" ;then
  echo "\$output6 Already exists" ;
  echo "Skipping Step5b for ${sample_id}" ;else 
  echo "Step 6 of Analysis Started for ${sample_id}"; 
  java -jar \${GATK} HaplotypeCaller -R \${ref_path}"\$reference" -I \$output5b -O \$output6 -ERC GVCF; 
  
    if test -f "\$output6"; then
      echo "\$output6 exists" ;
    else
      echo "Error!\$output6 not found";
      continue;
    fi
  fi
  
  
  module load zlib/1.2.8
  module load bcftools/1.4-6
  bcftools norm -m -any -Oz -o \${output6%.g.vcf.gz}"_temp.vcf.gz" \$output6
  bcftools view -e 'ALT="<NON_REF>"' -v snps,indels -Oz -o \${output6%.g.vcf.gz}"_mapped.bam.vcf.gz" \${output6%.g.vcf.gz}"_temp.vcf.gz"
  bcftools index -t \${output6%.g.vcf.gz}"_mapped.bam.vcf.gz"
  
  gunzip \${output6%.g.vcf.gz}"_mapped.bam.vcf.gz"
  
  
  echo "Step 6 of Analysis Completed for ${sample_id}"
  echo "END"
  
  
  
                                                          #STEP7
  module load gcc/12.3.0-gcc
  module load perl/5.34.0-intel
  
  
  output7=\${output6%.g.vcf.gz}".hg38_multianno.txt";
  echo "Checking for output file=\$output7 for step 7" ; 
  
  
  
  if test -f "\$output7"; then
  echo "\$output7 exists" ;
  echo "Skipping step 7 for ${sample_id}" ;
  echo "WES Analysis completed for ${sample_id}" 
  else
  echo "Step 7 of Analysis for ${sample_id} Started" 
  
  \$ref_path/annovar/table_annovar.pl \${output6%.g.vcf.gz}"_mapped.bam.vcf" \$ref_path/annovar/humandb -buildver hg38 -out \${output6%.g.vcf.gz} -arg '-splicing 15',,,,,,, -remove -protocol refGene,cytoBand,gnomad41_genome,gnomad41_exome,revel,dbscsnv11,intervar_20180118,clinvar_20250721 -operation g,r,f,f,f,f,f,f -nastring . -vcfinput;
  
  
  output7=\${output6%.g.vcf.gz}".hg38_multianno.txt";
    if test -f "\$output7"; then
     echo "\$output7 exists" ;
     echo "Step 7 of Analysis Completed for ${sample_id}" ;
    else
      echo "Error!\$output7 not found";
      continue
    fi
  fi
  echo "WES Anlaysis Completed for ${sample_id}" ;
  rm \${output6%.g.vcf.gz}"_temp.vcf.gz"
  
  echo "Date = ${params.date}" ;


  """
}

process geneFilter {
  input:
      path multianno_files
      path gene_panel_file
      path output_file
      
  output:
      path "${output_file}"
  
  publishDir (path: "${params.run_dir}/output/gene_filter", mode: 'move')
    
  script:
  """
  python3 ${params.gene_filter_script} \
    ${gene_panel_file} \
    "${params.keep_intronic}" \
    ${output_file} \
    ${multianno_files.join(' ')}
    
  """
}

workflow {
  
// This code assumes all fastq files are stored in directories that start with "Sample_" or some other prefix separated by an underscore
  Channel
    .fromPath("${params.sample_dir}/*", type: 'dir', checkIfExists: true) // emits each sample directory
    .map { dir ->
        def name = dir.getName() // extracts sample names
    
        // Keep only sample ID: ex. Sample_CML-012_345_678 becomes CML-012
        def sample_id = name.split('_')[1]

        def reads = file("${dir}/*.fastq.gz") // grabs the FASTQ files inside

        tuple(sample_id, name, reads)
    }
    .groupTuple()
    .map { sample_id, dir_names, reads_lists -> tuple(sample_id, dir_names, reads_lists.flatten())}
    .set { ch_samples }
  
  // Debugging option: Check the channel to make sure it's formatted correctly  
  //ch_samples.view()
  
  concatFastqs(ch_samples)
  
  fastqcArray(concatFastqs.out.combined_fqs)
  WESanalysis(concatFastqs.out.combined_fqs)

  fastqc_ch = fastqcArray.out.fastqc_zips.collect()
  multiQC(fastqc_ch)

  multianno_ch = WESanalysis.out.multianno.collect()
  // Defining channels to run geneFilter twice, once with each list of gene panels
  gene_panel_ch = Channel.from(file(params.gene_panel_file_1), file(params.gene_panel_file_2))
  filtered_output_ch = Channel.from(file(params.output_filename_1), file(params.output_filename_2))
  
  geneFilter(multianno_ch, gene_panel_ch, filtered_output_ch)

}