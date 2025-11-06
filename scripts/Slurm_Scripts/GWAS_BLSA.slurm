#!/bin/bash

##The below arguments are necessary for a job array in slurm.

#SBATCH --mail-user=alex.contreras@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=24GB
#SBATCH --array=1-36%36
#SBATCH --job-name=APOE4_Strat_BLSA
#SBATCH --output=/data/h_vmac/contra2/01_APOE4_Strat/logs/02_GWAS/BLSA_GWAS_%A_%a.out

module load PLINK/1.9b_5.2

#Setting up main directory.
dir=/data/h_vmac/contra2/01_APOE4_Strat
plots_output_dir=/data/h_vmac/contra2/01_APOE4_Strat/02_GWAS/plots
results_output_dir=/data/h_vmac/contra2/01_APOE4_Strat/02_GWAS/results
GWAMA_output_dir=/data/h_vmac/contra2/01_APOE4_Strat/02_GWAS/GWAMA_conversion

#Select the correct covariate directory
cov="1_cov_withComorbidities_50plus"

#Pulling cohort, race, phenotype, domain, comorbidities, age, analysis type
covar=$(cat $dir/01_CovPheno/$cov/BLSA*driver | awk -v line="$SLURM_ARRAY_TASK_ID" '{if (NR == line) print $1}') 
covar_no_ext="${covar%.txt}"

echo -e "\n$covar\n##Getting all necessary info from covariate file\n"

cohort=`echo $covar | cut -d"_" -f1` #ADNI
race=`echo $covar | cut -d"_" -f2` #NHB
pheno=`echo $covar | cut -d"_" -f3` #memslopes

comorbid=`echo $covar | cut -d"_" -f4` #WithComorbidities
age=`echo $covar | cut -d"_" -f5` #50plus
analysis=`echo $covar | cut -d"_" -f8 | cut -d"." -f1` #APOE4neg (gets rid of '.txt'')

#selecting the correct covariates
if [[ "$analysis" == "APOE4neg" || "$analysis" == "APOE4pos" ]]; then
  if [[ "$race" == "NHW" ]]; then
    echo "Analysis: $analysis, Race: $race"
    covars="Age, sex, NHW_PC1, NHW_PC2, NHW_PC3, NHW_PC4, NHW_PC5"
    datarace="NHW"
  elif [[ "$race" == "NHB" ]]; then
    echo "Analysis: $analysis, Race: $race"
    covars="Age, sex, NHB_PC1, NHB_PC2, NHB_PC3, NHB_PC4, NHB_PC5"
    datarace="AllRaces"
  fi
elif [[ "$analysis" == "ALL" ]]; then
  if [[ "$race" == "NHW" ]]; then
    echo "Analysis: $analysis, Race: $race"
    int_covars="Age, sex, apoe4pos, NHW_PC1, NHW_PC2, NHW_PC3, NHW_PC4, NHW_PC5"
    datarace="NHW"
  elif [[ "$race" == "NHB" ]]; then
    echo "Analysis: $analysis, Race: $race"
    int_covars="Age, sex, apoe4pos, NHB_PC1, NHB_PC2, NHB_PC3, NHB_PC4, NHB_PC5"
    datarace="AllRaces"
  fi
fi

################################################################################################################################################
# RUNNING GWAS
################################################################################################################################################

# Selecting the correct GWAS
if [[ "$analysis" == "APOE4neg" || "$analysis" == "APOE4pos" ]]; then
    echo "Performing GWAS analysis for ${analysis}"
     Code for the APOE4neg or APOE4pos case
    plink --bfile /data/h_vmac/hohmantj/Alex_GWAS_Stratified_APOE/BLSA_${datarace}_imputed_final_no_relateds \
    --keep $dir/01_CovPheno/$cov/$covar \
    --allow-no-sex \
    --pheno $dir/01_CovPheno/pheno/${cohort}_${race}_${pheno}_ALL.txt \
    --pheno-name ${pheno} \
    --covar $dir/01_CovPheno/$cov/$covar \
    --covar-name $covars \
    --no-const-covar \
    --linear hide-covar \
    --ci 0.95 \
    --out $results_output_dir/${covar_no_ext}
    
    #Plotting and GWAMA conversion for strat
    Rscript $dir/scripts/R_Scripts/generate_manhattan_qq_plots.R \
      $results_output_dir/${covar_no_ext}.assoc.linear \
      $plots_output_dir

    perl $dir/scripts/Perl_Scripts/plink2GWAMA_modified.pl \
      $results_output_dir/${covar_no_ext}.assoc.linear \
      /data/h_vmac/hohmantj/Alex_GWAS_Stratified_APOE/BLSA_${datarace}_imputed_final_no_relateds.frq \
      $GWAMA_output_dir/${covar_no_ext}.assoc.linear.GWAMA
elif [[ "$analysis" == "ALL" ]]; then
    echo "Performing GWAS analysis for ${analysis}"
    # Code here for the ALL case
    plink --bfile /data/h_vmac/hohmantj/Alex_GWAS_Stratified_APOE/BLSA_${datarace}_imputed_final_no_relateds \
    --keep $dir/01_CovPheno/$cov/$covar \
    --allow-no-sex \
    --pheno $dir/01_CovPheno/pheno/${cohort}_${race}_${pheno}_ALL.txt \
    --pheno-name ${pheno} \
    --covar $dir/01_CovPheno/$cov/$covar \
    --covar-name $int_covars \
    --no-const-covar \
    --linear interaction \
    --parameters 1-9,12 \
    --ci 0.95 \
    --out $results_output_dir/${covar_no_ext}_int

    #Subset the interaction term
    grep -E "ADDxapoe4pos|TEST" "$results_output_dir/${covar_no_ext}_int.assoc.linear" > "$results_output_dir/${covar_no_ext}_int_reduced.assoc.linear"
    
    Rscript $dir/scripts/R_Scripts/generate_manhattan_qq_plots.R \
      $results_output_dir/${covar_no_ext}_int_reduced.assoc.linear \
      $plots_output_dir
      
    perl $dir/scripts/Perl_Scripts/plink2GWAMA_modified.pl \
      $results_output_dir/${covar_no_ext}_int_reduced.assoc.linear \
      /data/h_vmac/hohmantj/Alex_GWAS_Stratified_APOE/BLSA_${datarace}_imputed_final_no_relateds.frq \
      $GWAMA_output_dir/${covar_no_ext}_int_reduced.assoc.linear.GWAMA
fi

