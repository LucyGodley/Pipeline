#PBS -N Step7
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=40gb
#PBS -o /scratch/path/logfiles/log.step7_new.out
#PBS -e /scratch/path/logfiles/log.step7_new.err
#PBS -d /scratch/path
module load gcc/4.9.4
module load perl/5.18.4
echo START
for input in `ls /scratch/path/*vcf`; do /gpfs/data/godley-lab/WES_analysis/annovar/table_annovar.pl $input /gpfs/data/godley-lab/WES_analysis/annovar/humandb/ -buildver hg38 -out ${input%.vcf} -arg '-splicing 15',,,,,,,, -remove -protocol refGene,cytoBand,exac03,gnomad_exome,gnomad_genome,kaviar_20150923,clinvar_20170905,avsnp147,dbnsfp30a -operation g,r,f,f,f,f,f,f,f -nastring . -vcfinput ; done
echo END
