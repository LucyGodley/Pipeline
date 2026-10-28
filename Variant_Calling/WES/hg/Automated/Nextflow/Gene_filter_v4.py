# Gene Filtration V4 based on Population Frequency
# Formatted to work with in-house nf-core pipeline
# Depth Cutoff 20X
# VAF Cutoff 30%
# Version 4.0.0
# Last edit: 7 Oct 2026 - Fixed incorrect column naming which resulted in incorrect AF values
  
import numpy as np
import pandas as pd
import glob
import os
from collections import defaultdict
import sys
  
vcf_columns_dict = {
  "v1": ['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','gnomAD_exome_ALL','gnomAD_exome_AFR','gnomAD_exome_AMR','gnomAD_exome_ASJ','gnomAD_exome_EAS','gnomAD_exome_FIN','gnomAD_exome_NFE','gnomAD_exome_OTH','gnomAD_exome_SAS','gnomAD_genome_ALL','gnomAD_genome_AFR','gnomAD_genome_AMR','gnomAD_genome_ASJ','gnomAD_genome_EAS','gnomAD_genome_FIN','gnomAD_genome_NFE','gnomAD_genome_OTH','Kaviar_AF','Kaviar_AC','Kaviar_AN','CLINSIG','CLNDBN','CLNACC','CLNDSDB','CLNDSDBID','avsnp147','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred','MutationAssessor_score','MutationAssessor_pred','FATHMM_score', 'FATHMM_pred','PROVEAN_score','PROVEAN_pred','VEST3_score','CADD_raw','CADD_phred','DANN_score','fathmm-MKL_coding_score','fathmm-MKL_coding_pred','MetaSVM_score','MetaSVM_pred','MetaLR_score','MetaLR_pred','integrated_fitCons_score','integrated_confidence_value','GERP++_RS','phyloP7way_vertebrate','phyloP20way_mammalian','phastCons7way_vertebrate','phastCons20way_mammalian','SiPhy_29way_logOdds','Otherinfo','#QUAL','DP','#CHROM','POS','#ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SeqInfo'],
  "v2": ['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','gnomAD_exome_ALL','AF_popmax','AF_male','AF_female','AF_raw','AF_afr','AF_sas','AF_amr','AF_eas','AF_nfe','AF_fin','AF_asj','AF_oth','non_topmed_AF_popmax','non_neuro_AF_popmax','non_cancer_AF_popmax','controls_AF_popmax','gnomAD_genome_ALL','gnomad312_AF_raw','gnomad312_AF_XX','gnomad312_AF_XY','gnomad312_AF_popmax','gnomad312_AF_faf95_popmax','gnomad312_AF_afr','gnomad312_AF_ami','gnomad312_AF_amr','gnomad312_AF_asj','gnomad312_AF_eas','gnomad312_AF_fin','gnomad312_AF_mid','gnomad312_AF_nfe','gnomad312_AF_oth','gnomad312_AF_sas','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','avsnp150','SIFT_score','SIFT_converted_rankscore','SIFT_pred','SIFT4G_score','SIFT4G_converted_rankscore','SIFT4G_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_rankscore','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_rankscore','Polyphen2_HVAR_pred','LRT_score','LRT_converted_rankscore','LRT_pred','MutationTaster_score','MutationTaster_converted_rankscore','MutationTaster_pred','MutationAssessor_score','MutationAssessor_rankscore','MutationAssessor_pred','FATHMM_score','FATHMM_converted_rankscore','FATHMM_pred','PROVEAN_score','PROVEAN_converted_rankscore','PROVEAN_pred','VEST4_score','VEST4_rankscore','MetaSVM_score','MetaSVM_rankscore','MetaSVM_pred','MetaLR_score','MetaLR_rankscore',\
  'MetaLR_pred','MetaRNN_score','MetaRNN_rankscore','MetaRNN_pred','M-CAP_score','M-CAP_rankscore','M-CAP_pred','REVEL_score','REVEL_rankscore','MutPred_score','MutPred_rankscore','MVP_score','MVP_rankscore','MPC_score','MPC_rankscore','PrimateAI_score','PrimateAI_rankscore','PrimateAI_pred','DEOGEN2_score','DEOGEN2_rankscore','DEOGEN2_pred','BayesDel_addAF_score','BayesDel_addAF_rankscore','BayesDel_addAF_pred','BayesDel_noAF_score','BayesDel_noAF_rankscore','BayesDel_noAF_pred','ClinPred_score','ClinPred_rankscore','ClinPred_pred','LIST-S2_score','LIST-S2_rankscore','LIST-S2_pred','Aloft_pred','Aloft_Confidence','CADD_raw','CADD_raw_rankscore','CADD_phred','DANN_score','DANN_rankscore','fathmm-MKL_coding_score','fathmm-MKL_coding_rankscore','fathmm-MKL_coding_pred','fathmm-XF_coding_score','fathmm-XF_coding_rankscore','fathmm-XF_coding_pred','Eigen-raw_coding','Eigen-raw_coding_rankscore','Eigen-PC-raw_coding','Eigen-PC-raw_coding_rankscore','GenoCanyon_score','GenoCanyon_rankscore','integrated_fitCons_score','integrated_fitCons_rankscore','integrated_confidence_value','LINSIGHT','LINSIGHT_rankscore','GERP++_NR','GERP++_RS','GERP++_RS_rankscore','phyloP100way_vertebrate','phyloP100way_vertebrate_rankscore','phyloP30way_mammalian','phyloP30way_mammalian_rankscore','phastCons100way_vertebrate','phastCons100way_vertebrate_rankscore','phastCons30way_mammalian','phastCons30way_mammalian_rankscore','SiPhy_29way_logOdds','SiPhy_29way_logOdds_rankscore','Interpro_domain','GTEx_V8_gene','GTEx_V8_tissue','Otherinfo','#QUAL','DP','#CHROM','POS','#ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SeqInfo'],
  "v3": ['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','gnomad41_exome_AF','gnomad41_exome_AF_raw','gnomad41_exome_AF_XX','gnomad41_exome_AF_XY','gnomad41_exome_AF_grpmax','gnomad41_exome_faf95','gnomad41_exome_faf99','gnomad41_exome_fafmax_faf95_max','gnomad41_exome_fafmax_faf99_max','gnomad41_exome_AF_afr','gnomad41_exome_AF_amr','gnomad41_exome_AF_asj','gnomad41_exome_AF_eas','gnomad41_exome_AF_fin','gnomad41_exome_AF_mid','gnomad41_exome_AF_nfe','gnomad41_exome_AF_remaining','gnomad41_exome_AF_sas','gnomad41_genome_AF','gnomad41_genome_AF_raw','gnomad41_genome_AF_XX','gnomad41_genome_AF_XY','gnomad41_genome_AF_grpmax','gnomad41_genome_faf95','gnomad41_genome_faf99','gnomad41_genome_fafmax_faf95_max','gnomad41_genome_fafmax_faf99_max','gnomad41_genome_AF_afr','gnomad41_genome_AF_ami','gnomad41_genome_AF_amr','gnomad41_genome_AF_asj','gnomad41_genome_AF_eas','gnomad41_genome_AF_fin','gnomad41_genome_AF_mid','gnomad41_genome_AF_nfe','gnomad41_genome_AF_remaining','gnomad41_genome_AF_sas','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','ONCDN','ONCDISDB','ONCREVSTAT','ONC','SCIDN','SCIDISDB','SCIREVSTAT','SCI','REVEL','InterVar_automated','PVS1','PS1','PS2','PS3','PS4','PM1','PM2','PM3','PM4','PM5','PM6','PP1','PP2','PP3','PP4','PP5','BA1','BS1','BS2','BS3','BS4','BP1','BP2','BP3','BP4','BP5','BP6','BP7','dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE','Otherinfo1','Otherinfo2','Otherinfo3',	'Otherinfo4',	'Otherinfo5','Otherinfo6','Otherinfo7','Otherinfo8','Otherinfo9','Otherinfo10','Otherinfo11','Otherinfo12','Otherinfo13']
}
    
keep_columns=['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','gnomad41_exome_AF','gnomad41_exome_AF_raw','gnomad41_exome_AF_XX','gnomad41_exome_AF_XY','gnomad41_exome_AF_grpmax','gnomad41_exome_faf95','gnomad41_exome_faf99','gnomad41_exome_fafmax_faf95_max','gnomad41_exome_fafmax_faf99_max','gnomad41_exome_AF_afr','gnomad41_exome_AF_amr','gnomad41_exome_AF_asj','gnomad41_exome_AF_eas','gnomad41_exome_AF_fin','gnomad41_exome_AF_mid','gnomad41_exome_AF_nfe','gnomad41_exome_AF_remaining','gnomad41_exome_AF_sas','gnomad41_genome_AF','gnomad41_genome_AF_raw','gnomad41_genome_AF_XX','gnomad41_genome_AF_XY','gnomad41_genome_AF_grpmax','gnomad41_genome_faf95','gnomad41_genome_faf99','gnomad41_genome_fafmax_faf95_max','gnomad41_genome_fafmax_faf99_max','gnomad41_genome_AF_afr','gnomad41_genome_AF_ami','gnomad41_genome_AF_amr','gnomad41_genome_AF_asj','gnomad41_genome_AF_eas','gnomad41_genome_AF_fin','gnomad41_genome_AF_mid','gnomad41_genome_AF_nfe','gnomad41_genome_AF_remaining','gnomad41_genome_AF_sas','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','ONCDN','ONCDISDB','ONCREVSTAT','ONC','SCIDN','SCIDISDB','SCIREVSTAT','SCI','REVEL','InterVar_automated','PVS1','PS1','PS2','PS3','PS4','PM1','PM2','PM3','PM4','PM5','PM6','PP1','PP2','PP3','PP4','PP5','BA1','BS1','BS2','BS3','BS4','BP1','BP2','BP3','BP4','BP5','BP6','BP7','dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE','Otherinfo','#QUAL','DP','#CHROM','POS','#ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SeqInfo']
  
  # Load gene panel
gene_panel_file = sys.argv[1]
with open(gene_panel_file, "r") as f:
    genes = f.read().split("\n")
  
high_pop_genes = ["CHEK2", "DDX41", "CFTR", "HFE"] #Added CFTR and HFE to catch AR carrier variants
keep_intronic = sys.argv[2]
keep_intronic = genes if keep_intronic == "Yes" else [""]
  
  # Filtering parameters
depth = 20  
vaf = 0.30
  
appended_data = []

output_filename = sys.argv[3]
  
df = pd.DataFrame(columns=['Parameters', 'Value'])
sample_list = []


def try_columns(file):
    # Read using the actual header in the file
    df = pd.read_csv(file, delimiter='\t', dtype=str, header=0)

    for key, columns in vcf_columns_dict.items():

        # Check whether this file matches any of the column sets defined above
        if set(columns).issubset(set(df.columns)):
            print(f"Using {key} column schema")

            # v3 uses Otherinfo1-Otherinfo13.
            # Rename them to the standard names used by v1/v2.
            if key == "v3":
                otherinfo_mapping = {
                    'Otherinfo1': 'Otherinfo',
                    'Otherinfo2': '#QUAL',
                    'Otherinfo3': 'DP',
                    'Otherinfo4': '#CHROM',
                    'Otherinfo5': 'POS',
                    'Otherinfo6': '#ID',
                    'Otherinfo7': 'REF',
                    'Otherinfo8': 'ALT',
                    'Otherinfo9': 'QUAL',
                    'Otherinfo10': 'FILTER',
                    'Otherinfo11': 'INFO',
                    'Otherinfo12': 'FORMAT',
                    'Otherinfo13': 'SeqInfo'
                }

                df = df.rename(columns=otherinfo_mapping)

            # Select the standardized columns
            return df[keep_columns]

    raise ValueError("File format does not match any known column schema.")

# ...existing code...
  
  
# Loop through each path
paths = sys.argv[4:]
for file in paths:
    print(file)
    filename = os.path.basename(file)
    sample = filename.split("_")[0]
    print(sample)
    sample_list.append(sample)
  

    # Read file using flexible columns
    unfiltered = try_columns(file)
    if '#QUAL' in unfiltered.columns:
        unfiltered.drop(columns=['#QUAL'], inplace=True, errors='ignore')
    if '#ID' in unfiltered.columns:
        unfiltered.drop(columns=['#ID'], inplace=True, errors='ignore')

    unfiltered['ID'] = sample
    if 'gnomad41_genome_AF' in unfiltered.columns:
             unfiltered['gnomad41_genome_AF'] = pd.to_numeric(
             unfiltered['gnomad41_genome_AF'].replace('.', 0),
             errors='coerce'
               ).fillna(0)
    else:
         unfiltered['gnomad41_genome_AF'] = 0

    pop_filtered_1 = unfiltered[unfiltered['gnomad41_genome_AF'] <= 0.01]
    pop_filtered_2 = unfiltered[unfiltered['gnomad41_genome_AF'] <= 0.05]

    gene_filtered_1 = pop_filtered_1[pop_filtered_1['Gene.refGene'].isin(genes)]
    gene_filtered_2 = pop_filtered_2[pop_filtered_2['Gene.refGene'].isin(high_pop_genes)]
    gene_filtered = pd.concat([gene_filtered_1, gene_filtered_2])

    gene_filtered = gene_filtered[gene_filtered["DP"] != "."]
    gene_filtered['AF'] = gene_filtered['SeqInfo'].str.split(':').str[1].str.split(',').str[-1]
    gene_filtered['AF'] = pd.to_numeric(gene_filtered['AF'], errors='coerce')
    gene_filtered['DP'] = pd.to_numeric(gene_filtered['DP'], errors='coerce')

    gene_filtered['VAF'] = gene_filtered['AF'] / gene_filtered['DP']
    gene_filtered = gene_filtered[gene_filtered['DP'] >= depth]
    gene_filtered = gene_filtered[gene_filtered['VAF'] >= vaf]

    # Remove duplicate entries
    gene_filtered = gene_filtered.drop_duplicates(keep='last')

    # Filter intronic variants if not needed
    mask = (gene_filtered['Func.refGene'] == "intronic") & (~gene_filtered['Gene.refGene'].isin(keep_intronic))
    intronic_filtered = gene_filtered[~mask].drop_duplicates(keep='last')

    appended_data.append(intronic_filtered)
        
# Parameter Summary
parameters = [
    ['Depth Cutoff', depth],
    ['VAF Cutoff', vaf],
    ['High Population frequency Genes(5%)', high_pop_genes],
    ['No of Genes', len(genes)],
    ['Gene List', genes],
    ['Number of Samples', len(sample_list)],
    ['Samples', sample_list]
]
for param in parameters:
    df.loc[len(df)] = param
    print(param)

# Save output
OUT = pd.concat(appended_data, sort=False)
with pd.ExcelWriter(output_filename) as writer:
    OUT.to_excel(writer, sheet_name='Filtered', index=False)
    df.to_excel(writer, sheet_name='Parameter', index=False)   
