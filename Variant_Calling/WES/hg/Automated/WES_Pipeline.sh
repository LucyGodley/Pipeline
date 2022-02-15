#PBS -N WES_Pipeline
#PBS -S /bin/bash
#PBS -l walltime=148:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=100gb
#PBS -o /path/logfiles/Output.out
#PBS -e /path/logfiles/Error.err
#PBS -d /path



logfile="/path/Logfile.txt"
mkdir Completed
path="/path/Completed"
echo "WES Analysis Started" > $logfile;

echo "Date = `date`" >> $logfile;

for input in `ls *_R1_001.fastq.gz`;do
SampleID=${input%_L*R1_001.fastq.gz};
lane=`ls ${SampleID}*R1*|wc -l`;

output7=`ls $SampleID*multianno.txt`
echo "Checking Final files=$output7 " >> $logfile;
if test -f "$output7"; then
  echo "$output7 exists! Moving the fastq files of $SampleID" >> $logfile;      
  mv `ls $SampleID*fastq.gz` $path;
  continue
else 
  echo "Running the WES pipeline" >> $logfile;
fi






                                                                          #STEP1
                                                             
module load java-jdk/1.8.0_92
module load picard/2.8.1


output1=${input%_R1_001.fastq.gz}"_ubam.bam"
echo "Checking output file=$output1 for step 1" >> $logfile;
if test -f "$output1"; then
  echo "$output1 Already exists" >> $logfile;
  echo "Skipping Step 1 for $input" >> $logfile;
else
echo "Step 1 of Analysis Started for $SampleID" >> $logfile
platform_unit=("`zless $input |head -n1 | awk -F ':' '{print $3 "." $4}'`");
java -Djava.io.tmpdir=${input%_R1_001.fastq.gz}"_001_tmp" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar FastqToSam FASTQ=$input FASTQ2=${input%_R1_001.fastq.gz}"_R2_001.fastq.gz" OUTPUT=${input%_R1_001.fastq.gz}"_ubam.bam" READ_GROUP_NAME=G SAMPLE_NAME=$SampleID LIBRARY_NAME=`(echo $input| awk -F '[_]' '{print $3}')` PLATFORM_UNIT=$platform_unit PLATFORM=illumina;

  if test -f "$output1"; then
    echo "$output1 exists" >> $logfile;
    echo "Step 1 of Analysis Completed for $SampleID" >> $logfile;
  else
     echo "Error! $output1 not found" >> $logfile; 
  fi
fi

module unload java-jdk/1.8.0_92
module unload picard/2.8.1    


                                                                       #STEP2
module load java-jdk/1.8.0_92
module load picard/2.8.1
output2=${output1%_ubam.bam}"_step2.bam";
echo "Checking for output file=$output2 for step 2" >> $logfile;

if test -f "$output2"; then
echo "$output2 Already exists" >> $logfile;
echo "Skipping Step 2 for $SampleID" >> $logfile;                                    else
echo "Step 2 of Analysis Started for $SampleID" >> $logfile;
java -Djava.io.tmpdir=${output1%_ubam.bam}"_tmp_step2" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MarkIlluminaAdapters I=$output1 O=${output1%_ubam.bam}"_step2.bam" M=${output1%_ubam.bam}"_step2_metrics.txt" ;
output2_metrics=${output1%_ubam.bam}"_step2_metrics.txt"
  if test -f "$output2"; then
     echo "$output2 exists" >> $logfile;
     echo "Step 2 of Analysis Completed for $SampleID" >> $logfile;
  else
    echo "Error! $output2 not found!" >> $logfile;
  fi
fi

module unload java-jdk/1.8.0_92
module unload picard/2.8.1 
 
                                                                     #STEP3A

module load java-jdk/1.8.0_92
module load picard/2.8.1

output3a=${output2%_step2.bam}"_3a_forBWA.fastq"
echo "Checking for output file=$output3a for step 3a" >> $logfile; 

if test -f "$output3a"; then 
echo "$output3a Already exists" >> $logfile;
echo "Skipping Step3a for $SampleID" >> $logfile;
else 
echo "Step 3 of Analysis Started for $SampleID" >> $logfile;
java -Djava.io.tmpdir=${output2%_step2.bam}"_tmp_step3a_BWA" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar SamToFastq I=$output2 FASTQ=${output2%_step2.bam}"_3a_forBWA.fastq" CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true;

  if test -f "$output3a"; then
     echo "$output3a exists" >> $logfile;
     echo "Step 3a of Analysis Completed for $SampleID" >> $logfile;
  else
     echo"Error! $output3a not found!">> $logfile;
  fi
fi

module unload java-jdk/1.8.0_92
module unload picard/2.8.1


                                                                   #STEP3B
module load gcc/6.2.0
module load java-jdk/1.8.0_92
module load picard/2.8.1
module load bwa/0.7.15


output3b=${output3a%_3a_forBWA.fastq}"_3b_bwa_mem.sam";
echo "Checking for output file=$output3b for step 3b" >> $logfile; 

if test -f "$output3b" ;then
echo "$output3b Already exists" >> $logfile;
echo "Skipping Step3b for $SampleID" >> $logfile;
else
echo "Step 3b of Analysis Started for $SampleID" >> $logfile;
/apps/software/gcc-6.2.0/bwa/0.7.15/bwa mem -M -t 7 -p /path/to/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz $output3a > ${output3a%_3a_forBWA.fastq}"_3b_bwa_mem.sam" ;
  if test -f "$output3b"; then
    echo "$output3b exists" >> $logfile;
    echo "Step 3b of Analysis Completed for $SampleID" >> $logfile;
  else
    echo "Error!$output3b not found!">> $logfile;
    
  fi
fi


module unload gcc/6.2.0
module unload java-jdk/1.8.0_92
module unload picard/2.8.1
module unload bwa/0.7.15


                                                                 #STEP3C    
module load java-jdk/1.8.0_92
module load picard/2.8.1
output3c=${output3b%_3b_bwa_mem.sam}"_mapped.bam"
echo "Checking for output file=$output3c for step 3c" >> $logfile;

if test -f "$output3c" ;then
echo "$output3c Already exists" >> $logfile;
echo "Skipping Step3c for $SampleID" >> $logfile;
else
echo "Step 3c of Analysis Started for $SampleID" >> $logfile; 
java -Djava.io.tmpdir=${output3b%_3b_bwa_mem.sam}"_tmp_step3c_BWA" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MergeBamAlignment R=/path/to/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz ALIGNED_BAM=$output3b UNMAPPED_BAM=${output3b%_3b_bwa_mem.sam}"_ubam.bam" OUTPUT=${output3b%_3b_bwa_mem.sam}"_mapped.bam" CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS;


  if test -f "$output3c"; then
    echo "$output3c exists" >> $logfile;
    echo "Step3c of Analysis Completed for $SampleID" >> $logfile;
  else
    echo "Error! $output3c not found!" >> $logfile;
  fi
fi


module unload java-jdk/1.8.0_92
module unload picard/2.8.1
 
                                                                #STEP4 
module load gcc/6.2.0
module load samtools/1.9


output4=`ls $SampleID*_step4.bam`;
echo "Checking for output file=$output4 for step 4" >> $logfile ;
if test -f "$output4" ;then
echo "$output4 Already exists" >> $logfile;
echo "Skipping Step4 for $SampleID" >> $logfile;
else
echo "Step 4 of Analysis Started for $SampleID" >> $logfile;
case $lane in
   1) echo "$SampleID has 1 Lanes" >> $logfile;
files=($(ls -d $SampleID*_mapped.bam));
echo Number of Mapped Bam files= "${#files[@]}">> $logfile;
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MarkDuplicates INPUT=$output3c OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true; 
output4=${output3c%$detach}"_step4.bam";   
;; 
   2) echo "$SampleID has 2 Lanes"  >> $logfile;
   if [ `ls ${SampleID}*mapped.bam|wc -l` -eq $lane ];then
files=($(ls -d $SampleID*_mapped.bam));
echo Number of Mapped Bam files= "${#files[@]}">> $logfile;
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MarkDuplicates INPUT=${files[0]} INPUT=${files[1]} OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true;
output4=${output3c%$detach}"_step4.bam"  
 else
 continue 
 fi
 ;;
3) echo "$SampleID has 3 Lanes" >> $logfile
if [ `ls ${SampleID}*mapped.bam|wc -l` -eq $lane ];then
files=($(ls -d $SampleID*_mapped.bam))
echo Number of Mapped Bam files= "${#files[@]}">> $logfile;
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MarkDuplicates INPUT=${files[0]} INPUT=${files[1]} INPUT=${files[2]} OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true;
output4=${output3c%$detach}"_step4.bam";
else
continue;
fi   
;;
4) echo "$SampleID has 4 Lanes" >> $logfile
if [ `ls ${SampleID}*mapped.bam|wc -l` -eq $lane ];then
files=($(ls -d $SampleID*mapped.bam))
echo Number of Mapped Bam files= "${#files[@]}">> $logfile;
detach=`(echo $output3c| awk -F '[_]' '{print "_" $3 "_mapped.bam"}')`
java -Djava.io/.tmpdir=${output3c%_mapped.bam}"_tmp_step4" -jar /apps/software/java-jdk-1.8.0_92/picard/2.8.1/picard.jar MarkDuplicates INPUT=${files[0]} INPUT=${files[1]} INPUT=${files[2]} INPUT=${files[3]} OUTPUT=${output3c%$detach}"_step4.bam" METRICS_FILE=${output3c%$detach}"_step4_metrics.txt" OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true; 
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
    echo "$output4 exists" >> $logfile;
    echo "Step4 of Analysis Completed for $SampleID" >> $logfile;
  else
    echo "Error! $output4 not found">> $logfile;
  fi
fi

module unload java-jdk/1.8.0_92
module unload picard/2.8.1 

                                                            #STEP5a


module load java-jdk/1.8.0_92
module load picard/2.8.1
output5a=${output4%_step4.bam}"_step5a_BaseRecalibrator.table";
echo "Checking for output file=$output5a for step 5a" >> $logfile;

if test -f "$output5a" ;then
echo "$output5a Already exists" >> $logfile;
echo "Skipping Step5a for $SampleID" >> $logfile;
else
echo "Step 5a of Analysis for $SampleID Started" >> $logfile;
java -jar /apps/software/java-jdk-1.8.0_92/gatk/3.7/GenomeAnalysisTK.jar -T BaseRecalibrator -R /path/to/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa -I $output4 -knownSites /path/to/dbsnp/All_20170710_sort.vcf -o ${output4%_step4.bam}"_step5a_BaseRecalibrator.table" -U ALLOW_SEQ_DICT_INCOMPATIBILITY; 

  if test -f "$output5a"; then
    echo "$output5a exists" >> $logfile;
    echo "Step 5a of Analysis Completed" >> $logfile;
  else
    echo "Error!$output5a not found!">> $logfile;
  fi
fi


module unload java-jdk/1.8.0_92
module unload picard/2.8.1


                                                         #STEP5b

module load java-jdk/1.8.0_92
module load picard/2.8.1

output5b=${output4%_step4.bam}"_step5b.bam";
echo "Checking for output file=$output5b of step 5b" >> $logfile; 

if test -f "$output5b" ;then
echo "$output5b Already exists" >> $logfile;
echo "Skipping Step5b for $SampleID" >> $logfile;                                                                           else 
echo "Step 5b of Analysis Started" >> $logfile;
java -jar /apps/software/java-jdk-1.8.0_92/gatk/3.7/GenomeAnalysisTK.jar -T PrintReads -R /path/to/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa -I $output4 -BQSR ${output4%_step4.bam}"_step5a_BaseRecalibrator.table" -o ${output4%_step4.bam}"_step5b.bam" ;

  if test -f "$output5b"; then
    echo "$output5b exists" >> $logfile;
    echo "Step 5b of Analysis Completed for $SampleID" >> $logfile;
  else
    echo "Error! $output5b not found! ">> $logfile;
  fi
fi
 
module unload java-jdk/1.8.0_92
module unload picard/2.8.1


                                                          #STEP6

module load java-jdk/1.8.0_92
module load gatk/3.7
output6=${output5b%_step5b.bam}".vcf"
echo "Checking for output file=$output6 for step 6" >> $logfile; 

if test -f "$output6" ;then
echo "$output6 Already exists" >> $logfile;
echo "Skipping Step5b for $SampleID" >> $logfile;                                                                           else 
echo "Step 6 of Analysis Started for $SampleID" >>$logfile
java -jar /apps/software/java-jdk-1.8.0_92/gatk/3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R /path/to/reference_human_38/Homo_sapiens.GRCh38.dna.toplevel.fa -I $output5b --genotyping_mode DISCOVERY -o ${output5b%_step5b.bam}".vcf"; 

  if test -f "$output6"; then
    echo "$output6 exists" >> $logfile;
    echo "Step 6 of Analysis Completed for $SampleID" >> $logfile;
  else
    echo "Error!$output6 not found">> $logfile;
  fi
fi

module unload java-jdk/1.8.0_92
module unload gatk/3.7 

                                                        #STEP7
module load gcc/4.9.4
module load perl/5.18.4


output7=`ls $SampleID*multianno.txt`;
echo "Checking for output file=$output7 for step 7" >> $logfile; 



if test -f "$output7"; then
echo "$output7 exists" >> $logfile;
echo "Skipping step 7 for $SampleID" >> $logfile;
echo "WES Analysis completed for $SampleID" >> $logfile
else
echo "Step 7 of Analysis for $SampleID Started" >> $logfile
/gpfs/data/godley-lab/WES_analysis/annovar/table_annovar.pl $output6 /path/to/annovar/humandb/ -buildver hg38 -out ${output6%.vcf} -arg '-splicing 15',,,,,,,, -remove -protocol refGene,cytoBand,exac03,gnomad_exome,gnomad_genome,kaviar_20150923,clinvar_20170905,avsnp147,dbnsfp30a -operation g,r,f,f,f,f,f,f,f -nastring . -vcfinput ; 

output7=`ls $SampleID*multianno.txt`;
  if test -f "$output7"; then
    echo "$output7 exists" >> $logfile;
    echo "Step 7 of Analysis Completed for $SampleID" >> $logfile;
    mv `ls $SampleID*fastq.gz` $path 
  else
    echo "Error!$output7 not found">> $logfile;
  fi
fi
echo "WES Anlaysis Completed for $SampleID" >> $logfile;
done
echo "Date = `date`" >> $logfile;
