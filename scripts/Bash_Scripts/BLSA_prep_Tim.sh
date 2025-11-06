#!/bin/bash

# Create directory
mkdir -p /data/h_vmac/hohmantj/Alex_GWAS_Stratified_APOE

module load PLINK/1.9b_5.2
# Run PLINK for BLSA_NHW_imputed_final
plink --bfile /data/h_vmac/hohmantj/BLSA/GWAS/Imputed/TOPMed/Cleaned/BLSA_NHW_imputed_final \
      --keep /data/h_vmac/contra2/CODE_CHECK/0_Genetic_Files/Data/NHW/BLSA_NHW_imputed_final_no_relateds.txt \
      --make-bed \
      --out /data/h_vmac/hohmantj/Alex_GWAS_Stratified_APOE/BLSA_NHW_imputed_final_no_relateds \
      --freq

# Run PLINK for BLSA_AllRaces_imputed_final
plink --bfile /data/h_vmac/hohmantj/BLSA/GWAS/Imputed/TOPMed/Cleaned/BLSA_AllRaces_imputed_final \
      --keep /data/h_vmac/contra2/CODE_CHECK/0_Genetic_Files/Data/AllRaces/BLSA_AllRaces_imputed_final_no_relateds.txt \
      --make-bed \
      --out /data/h_vmac/hohmantj/Alex_GWAS_Stratified_APOE/BLSA_AllRaces_imputed_final_no_relateds \
      --freq

