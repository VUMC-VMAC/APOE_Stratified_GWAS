#!/bin/bash
#SBATCH--mail-user=alex.contreras@vanderbilt.edu
#SBATCH--mail-type=ALL
#SBATCH--ntasks=1
#SBATCH--cpus-per-task=1
#SBATCH--time=80:00:00
#SBATCH--mem=24GB
#SBATCH--job-name=MAGMA_memslopes
#SBATCH--output=/data/h_vmac/contra2/01_APOE4_Strat/logs/05_MAGMA/MAGMA_memslopes_%A_%a.out

#Directories
magma_dir=/data/h_vmac/contra2/01_APOE4_Strat/05_MAGMA
input_dir=/data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/results
R_scripts_dir=/data/h_vmac/contra2/01_APOE4_Strat/scripts/R_Scripts
genetic_dir=/data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files

datarace=AllRaces
domain=memslopes

#file could be :  /data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/results/MetaCrossAnc_CrossAnc_memslopes_ALL.in_2of2Ancestry.out

# Loop through each file containing "_memslopes_" in the input directory
for file in "$input_dir"/*_${domain}_*_2of2Ancestry.out; do
    filename=$(basename "$file")
    basename="${filename%.*}"

    echo -e "\nProcessing file: $filename"

    # Running Gene Test
    echo -e "\nGene Test: $filename"
    magma --bfile $genetic_dir/${datarace}/${datarace}_combined \
          --pval "$file" use=rs_number,p-value ncol=n_samples \
          --gene-annot $magma_dir/annotation_files/${datarace}.genes.annot \
          --out $magma_dir/gene_results/${basename}_gene

    # Pathway test
    echo -e "\nPathway Test: $filename"
    magma --gene-results "$magma_dir/gene_results/${basename}_gene.genes.raw" \
          --set-annot /data/h_vmac/Programs/magma/custom_pathways.txt \
          --out $magma_dir/pathway_results/${basename}_pathway

    # FDR correction
    echo -e "\nFDR Correction: $filename"
    Rscript $R_scripts_dir/MAGMA_FDR_correction_and_threshold.R \
            $magma_dir/gene_results/${basename}_gene.genes.out \
            $magma_dir/pathway_results/${basename}_pathway.gsa.out

    # Manhattan Plots
    echo -e "\nManhattan Plots: $filename"
    Rscript $R_scripts_dir/MAGMA_gene_manhattan.R \
            $magma_dir/gene_results/${basename}_gene.genes.out
done