#!/bin/bash
#Do this in the directory with all the GWAS outputs (after GWAMA conversion) to make .in files
pushd /data/h_vmac/contra2/01_APOE4_Strat/02_GWAS/GWAMA_conversion

# List of all phenotypes
phenotypes=("NHW_MEM" "NHW_memslopes" "NHW_EXF" "NHW_exfslopes" "NHW_LAN" "NHW_lanslopes" "NHB_MEM" "NHB_memslopes" "NHB_EXF" "NHB_exfslopes" "NHB_LAN" "NHB_lanslopes")
# Loop through each phenotype
for phenotype in "${phenotypes[@]}"; do
    # Use ls and grep to filter filenames for APOE4neg and store in a list file
    ls | grep "${phenotype}.*APOE4neg.*GWAMA" > "${phenotype}_APOE4neg.in"
    # Use ls and grep to filter filenames for APOE4pos and store in a list file
    ls | grep "${phenotype}.*APOE4pos.*GWAMA" > "${phenotype}_APOE4pos.in"
    # Use ls and grep to filter filenames for ALL and store in a list file
    ls | grep "${phenotype}.*ALL.*GWAMA" > "${phenotype}_ALL.in"
done

ls *in > fileset_GWAMA_driver.txt

popd