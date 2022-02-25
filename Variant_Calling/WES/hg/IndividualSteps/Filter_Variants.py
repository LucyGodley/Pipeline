import numpy as np
import pandas as pd
import glob
import os
from collections import defaultdict

vcf_columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
       'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
       'cytoBand', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN',
       'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'gnomAD_exome_ALL',
       'gnomAD_exome_AFR', 'gnomAD_exome_AMR', 'gnomAD_exome_ASJ',
       'gnomAD_exome_EAS', 'gnomAD_exome_FIN', 'gnomAD_exome_NFE',
       'gnomAD_exome_OTH', 'gnomAD_exome_SAS', 'gnomAD_genome_ALL',
       'gnomAD_genome_AFR', 'gnomAD_genome_AMR', 'gnomAD_genome_ASJ',
       'gnomAD_genome_EAS', 'gnomAD_genome_FIN', 'gnomAD_genome_NFE',
       'gnomAD_genome_OTH', 'Kaviar_AF', 'Kaviar_AC', 'Kaviar_AN', 'CLINSIG',
       'CLNDBN', 'CLNACC', 'CLNDSDB', 'CLNDSDBID', 'avsnp147', 'SIFT_score',
       'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred',
       'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred',
       'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score',
       'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score',
       'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score',
       'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score',
       'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred',
       'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS',
       'phyloP7way_vertebrate', 'phyloP20way_mammalian',
       'phastCons7way_vertebrate', 'phastCons20way_mammalian',
       'SiPhy_29way_logOdds', 'Otherinfo', 'Unnamed: 80', 'Unnamed: 81',
       'Unnamed: 82', 'Unnamed: 83', 'Unnamed: 84', 'Unnamed: 85',
       'Unnamed: 86', 'Unnamed: 87', 'Unnamed: 88', 'Unnamed: 89',
       'Unnamed: 90', 'Unnamed: 91']

paths = []
gene_panel_file=input("Enter the file name of Gene Panel:\n")


with open(gene_panel_file,"r") as f:
 genes = f.read().split("\n")

keep_intronic=[]
keep_intronic= genes
print(keep_intronic)



appended_data = []


#get number of directories and paths
folders = input("How many folders do you have filtered files in? \n")
for i in range(int(folders)):
    path = input("Enter path to the folder %s the annotated_txt files:\n" %(i+1,))
    paths.append(path)

output_filename=input("Enter output filename with .xlsx at the end:\n")
print(output_filename)

print(genes)
#get number of genes and gene name
#gene_number = input("How many genes are you filtering?\n")
#for i in range(int(gene_number)):
 #   gene = input("Enter the name of gene  %s:\n" %(i+1,))
 #  genes.append(gene)
 #   intronic = input("Do you want to keep intronic variants for this gene? Enter y or n:\n")
  #  if intronic == "y":
     #   

i=0
for path in paths:
    for filename in glob.glob(os.path.join(path, '*multianno.txt')):
        i+=1
        sample = filename.split("/")[-1]
        #read in txt as a dataframe
        unfiltered = pd.read_csv(filename, delimiter = '\t',dtype = str, names = vcf_columns, lineterminator = "\n", header = 0)
        unfiltered["ID"] = sample

        #replace null data
        unfiltered["gnomAD_genome_ALL"] = unfiltered["gnomAD_genome_ALL"].replace('.', 0)
        #convert population frequencies to integer values in order to filter
        unfiltered["gnomAD_genome_ALL"] = unfiltered["gnomAD_genome_ALL"].apply(pd.to_numeric)


        #select rows from unfiltered where the string based on mutation type and population frequency
        #NOTE: If you dont want a population filter, or want a different population, change this line - can add a # to comment it out
        pop_filtered = unfiltered.loc[(unfiltered['gnomAD_genome_ALL'] <= 0.01)]

        # filter by genes of interest
        gene_filtered = pop_filtered.loc[pop_filtered['Gene.refGene'].isin(genes)]

        #keep/filter out intronic variants
        mask = (gene_filtered['Func.refGene'] == "intronic") & (~gene_filtered['Gene.refGene'].isin(keep_intronic))
        intronic_filtered = gene_filtered[~mask]
    
        appended_data.append(intronic_filtered)


#print(intronic_filtered)

print(i)
#OUT = pd.concat(appended_data, sort = False)
OUT = pd.concat(appended_data, sort = False)

#NOTE: change the path in output file for your own computer
output_file = [path+"/"+output_filename]
writer = pd.ExcelWriter(output_file[0])
OUT.to_excel(writer,'Filtered')
writer.save()
