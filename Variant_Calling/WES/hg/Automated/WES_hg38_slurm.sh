#!/bin/bash
#SBATCH -A p32126
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 96:00:00
#SBATCH --mem=50G
#SBATCH -J NCH_b7
#SBATCH -o NCH_b7_%j.out
#SBATCH -e NCH_b7_%j.err

###################################################

#UKBIOBANK WES PIPELINE
#DATE=10/03/2023
#AUTHOR=ASHWIN
#Removed additional log file
###################################################

ref_path="/projects/p32126/WES_analysis"
reference="/reference_hg_38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
dbsnp="/dbsnp/hg38/ALL_20141222.dbSNP142_human_GRCh38_sort.snps.vcf"

mkdir Completed

echo START
echo "WES Analysis Started" ;

echo "Working directory=`pwd`";
echo "Date = `date`" ;



for input in `ls *_R1_001.fastq.gz`;do
SampleID=$(echo "$input" |cut -d"_" -f1);
lane=`ls ${SampleID}*R1*|wc -l`;

echo $SampleID

output7=`ls $SampleID*multianno.txt`
echo "Checking Final files=$output7 " ;
if test -f "$output7"; then
  echo "$output7 exists! Moving the fastq files of $SampleID" ;      
  mv `ls $SampleID*fastq.gz` $path;
  continue
else 
  echo "Running the WES pipeline" ;
fi






                                                                          #STEP1
                                                             
module load picard/2.6.0
PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"

output1=${input%_R1_001.fastq.gz}"_ubam.bam"
echo "Checking output file=$output1 for step 1" ;
if test -f "$output1"; then
  echo "$output1 Already exists" ;
  echo "Skipping Step 1 for $input" ;
else
echo "Step 1 of Analysis Started for $SampleID" 
platform_unit=("`zless $input |head -n1 | awk -F ':' '{print $3 "." $4}'`");
java -Djava.io.tmpdir=${input%_R1_001.fastq.gz}"_001_tmp" -jar ${PICARD} FastqToSam FASTQ=$input FASTQ2=${input%_R1_001.fastq.gz}"_R2_001.fastq.gz" OUTPUT=${input%_R1_001.fastq.gz}"_ubam.bam" READ_GROUP_NAME=DR SAMPLE_NAME=$SampleID LIBRARY_NAME=`(echo $input| awk -F '[_]' '{print $(NF-2)}')` PLATFORM_UNIT=$platform_unit PLATFORM=ILLUMINA;

  if test -f "$output1"; then
    echo "$output1 exists" ;
    echo "Step 1 of Analysis Completed for $SampleID" ;
  else
     echo "Error! $output1 not found" ;
     continue; 
  fi
fi

  


                                                                       #STEP2

module load picard/2.6.0
PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"


output2=${output1%_ubam.bam}"_step2.bam";
echo "Checking for output file=$output2 for step 2" ;

if test -f "$output2"; then
echo "$output2 Already exists" ;
echo "Skipping Step 2 for $SampleID" ;                                    else
echo "Step 2 of Analysis Started for $SampleID" ;
java -Djava.io.tmpdir=${output1%_ubam.bam}"_tmp_step2" -jar ${PICARD} MarkIlluminaAdapters I=$output1 O=${output1%_ubam.bam}"_step2.bam" M=${output1%_ubam.bam}"_step2_metrics.txt" ;
output2_metrics=${output1%_ubam.bam}"_step2_metrics.txt"
  if test -f "$output2"; then
     echo "$output2 exists" ;
     echo "Step 2 of Analysis Completed for $SampleID" ;
  else
    echo "Error! $output2 not found!" ;
    continue;
  fi
fi


 
                                                                     #STEP3A


module load picard/2.6.0
PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"

output3a=${output2%_step2.bam}"_3a_forBWA.fastq"
echo "Checking for output file=$output3a for step 3a" ; 

if test -f "$output3a"; then 
echo "$output3a Already exists" ;
echo "Skipping Step3a for $SampleID" ;
else 
echo "Step 3 of Analysis Started for $SampleID" ;
java -Djava.io.tmpdir=${output2%_step2.bam}"_tmp_step3a_BWA" -jar ${PICARD} SamToFastq I=$output2 FASTQ=${output2%_step2.bam}"_3a_forBWA.fastq" CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true;

  if test -f "$output3a"; then
     echo "$output3a exists" ;
     echo "Step 3a of Analysis Completed for $SampleID" ;
  else
     echo"Error! $output3a not found!";
     continue;
  fi
fi




                                                                   #STEP3B

module load bwa/0.7.17
BWA="/software/bwa/gcc/0.7.17/bwa"

output3b=${output3a%_3a_forBWA.fastq}"_3b_bwa_mem.sam";
echo "Checking for output file=$output3b for step 3b" ; 

if test -f "$output3b" ;then
echo "$output3b Already exists" ;
echo "Skipping Step3b for $SampleID" ;
else
echo "Step 3b of Analysis Started for $SampleID" ;
${BWA} mem -M -t 7 -p ${ref_path}"$reference" $output3a > ${output3a%_3a_forBWA.fastq}"_3b_bwa_mem.sam" ;
  if test -f "$output3b"; then
    echo "$output3b exists" ;
    echo "Step 3b of Analysis Completed for $SampleID" ;
  else
    echo "Error!$output3b not found!";
    continue;
  fi
fi




                                                                 #STEP3C    

module load picard/2.6.0
PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"



output3c=${output3b%_3b_bwa_mem.sam}"_mapped.bam"
echo "Checking for output file=$output3c for step 3c" ;

if test -f "$output3c" ;then
echo "$output3c Already exists" ;
echo "Skipping Step3c for $SampleID" ;
else
echo "Step 3c of Analysis Started for $SampleID" ; 
java -Djava.io.tmpdir=${output3b%_3b_bwa_mem.sam}"_tmp_step3c_BWA" -jar ${PICARD} MergeBamAlignment R=${ref_path}"$reference" ALIGNED_BAM=$output3b UNMAPPED_BAM=${output3b%_3b_bwa_mem.sam}"_ubam.bam" OUTPUT=${output3b%_3b_bwa_mem.sam}"_mapped.bam" CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS;


  if test -f "$output3c"; then
    echo "$output3c exists" ;
    echo "Step3c of Analysis Completed for $SampleID" ;
  else
    echo "Error! $output3c not found!" ;
    continue;
  fi
fi



                                                                #STEP4 
module load picard/2.6.0
PICARD="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar"

output4=`ls $SampleID*_step4.bam`;
echo "Checking for output file=$output4 for step 4"  ;
if test -f "$output4" ;then
echo "$output4 Already exists" ;
echo "Skipping Step4 for $SampleID" ;
rm $SampleID*sam ;                                                                                                      rm $SampleID*fastq ;                                                                                                    rm $SampleID*ubam ;                                                                                                     rm $SampleID*step2* ;
else
echo "Step 4 of Analysis Started for $SampleID" ;
case $lane in
   1) echo "$SampleID has 1 Lanes" ;
files=($(ls -d $SampleID*_mapped.bam));
echo Number of Mapped Bam files= "${#files[@]}";
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar ${PICARD} MarkDuplicates INPUT=$output3c OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true; 
output4=${output3c%$detach}"_step4.bam";   
;; 
   2) echo "$SampleID has 2 Lanes"  ;
   if [ `ls ${SampleID}*mapped.bam|wc -l` -eq $lane ];then
files=($(ls -d $SampleID*_mapped.bam));
echo Number of Mapped Bam files= "${#files[@]}";
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar ${PICARD} MarkDuplicates INPUT=${files[0]} INPUT=${files[1]} OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true;
output4=${output3c%$detach}"_step4.bam"  
 else
 continue 
 fi
 ;;
3) echo "$SampleID has 3 Lanes" 
if [ `ls ${SampleID}*mapped.bam|wc -l` -eq $lane ];then
files=($(ls -d $SampleID*_mapped.bam))
echo Number of Mapped Bam files= "${#files[@]}";
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar ${PICARD} MarkDuplicates INPUT=${files[0]} INPUT=${files[1]} INPUT=${files[2]} OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true;
output4=${output3c%$detach}"_step4.bam";
else
continue;
fi   
;;
4) echo "$SampleID has 4 Lanes" 
if [ `ls ${SampleID}*mapped.bam|wc -l` -eq $lane ];then
files=($(ls -d $SampleID*mapped.bam))
echo Number of Mapped Bam files= "${#files[@]}";
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar ${PICARD} MarkDuplicates INPUT=${files[0]} INPUT=${files[1]} INPUT=${files[2]} INPUT=${files[3]} OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true; 
output4=${output3c%$detach}"_step4.bam";
else
continue
fi
 ;;
*)
 echo "Flop"
;;
esac

  if test -f "$output4"; then
    echo "$output4 exists" ;
    echo "Step4 of Analysis Completed for $SampleID" ;
    rm $SampleID*sam ;
    rm $SampleID*fastq ;
    rm $SampleID*ubam ;
    rm $SampleID*step2* ;
  else
    echo "Error! $output4 not found";
    continue;
  fi
fi


                                                            #STEP5a

module load java/jdk-17.0.2+8
module load gatk/4.4.0.0
GATK="/software/gatk/4.4.0.0/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"

output5a=${output4%_step4.bam}"_step5a_BaseRecalibrator.table";
echo "Checking for output file=$output5a for step 5a" ;

if test -f "$output5a" ;then
echo "$output5a Already exists" ;
echo "Skipping Step5a for $SampleID" ;
else
echo "Step 5a of Analysis for $SampleID Started" ;
java -jar ${GATK} BaseRecalibrator -R ${ref_path}"$reference" -I $output4 --known-sites ${ref_path}"$dbsnp" -O $output5a
  if test -f "$output5a"; then
    echo "$output5a exists" ;
    echo "Step 5a of Analysis Completed" ;
  else
    echo "Error!$output5a not found!";
    continue
  fi
fi





                                                         #STEP5b

module load java/jdk-17.0.2+8
module load gatk/4.4.0.0
GATK="/software/gatk/4.4.0.0/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"

output5b=${output4%_step4.bam}"_step5b.bam";
echo "Checking for output file=$output5b of step 5b" ; 

if test -f "$output5b" ;then
echo "$output5b Already exists" ;
echo "Skipping Step5b for $SampleID" ;                                                                           else 
echo "Step 5b of Analysis Started" ;
java -jar ${GATK} ApplyBQSR -R ${ref_path}"$reference" -I $output4 -bqsr $output5a -O ${output4%_step4.bam}"_step5b.bam" ;

  if test -f "$output5b"; then
    echo "$output5b exists" ;
    echo "Step 5b of Analysis Completed for $SampleID" ;
  else
   echo "Error! $output5b not found! ";
   continue;
  fi
fi
 




                                                        #STEP6

module load java/jdk-17.0.2+8
module load gatk/4.4.0.0
GATK="/software/gatk/4.4.0.0/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"

output6=${output5b%_step5b.bam}".vcf"
echo "Checking for output file=$output6 for step 6" ; 

if test -f "$output6" ;then
echo "$output6 Already exists" ;
echo "Skipping Step5b for $SampleID" ;                                                                           else 
echo "Step 6 of Analysis Started for $SampleID" 
java -jar ${GATK} HaplotypeCaller -R ${ref_path}"$reference" -I $output5b -O $output6; 

  if test -f "$output6"; then
    echo "$output6 exists" ;
    echo "Step 6 of Analysis Completed for $SampleID" ;
  else
    echo "Error!$output6 not found";
    continue;
  fi
fi

#gzip ${SampleID}*vcf*






                                                        #STEP7
module load gcc/12.3.0-gcc
module load perl/5.34.0-intel


output7=`ls $SampleID*multianno.txt`;
echo "Checking for output file=$output7 for step 7" ; 



if test -f "$output7"; then
echo "$output7 exists" ;
echo "Skipping step 7 for $SampleID" ;
echo "WES Analysis completed for $SampleID" 
else
echo "Step 7 of Analysis for $SampleID Started" 
$ref_path/annovar/table_annovar.pl $output6 $ref_path/annovar/humandb -buildver hg38 -out ${output6%.vcf} -arg '-splicing 15',,,,,,, -remove -protocol refGene,cytoBand,exac03,gnomad211_exome,gnomad312_genome,clinvar_20221231,avsnp150,dbnsfp42a -operation g,r,f,f,f,f,f,f -nastring . -vcfinput ; 

output7=`ls $SampleID*multianno.txt`;
  if test -f "$output7"; then
   echo "$output7 exists" ;
   echo "Step 7 of Analysis Completed for $SampleID" ;
    mv `ls $SampleID*fastq.gz` $path
  else
    echo "Error!$output7 not found";
    continue
  fi
fi
mv ${SampleID}*fastq.gz Completed/
echo "WES Anlaysis Completed for $SampleID" ;
echo END
done
echo "Date = `date`" ;

