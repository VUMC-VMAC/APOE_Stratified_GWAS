#!/bin/bash

#Make sure you are in the GNOVA environment before running this
pushd /data/h_vmac/contra2/01_APOE4_Strat/03_GWAMA/results

# Loop through each .out file in the NHW GWAMA sumstats directory 
for file in ./*NHW*cohorts.out; do 
    # Extract the filename 
    file_name=$(basename "$file" .out)
    echo -e "\nProcessing file: $file" 
    python2 /data/h_vmac/contra2/01_APOE4_Strat/scripts/Python_Scripts/GNOVA-2.0/munge_sumstats.py \
        --signed-sumstats beta,0 \
        --out "/data/h_vmac/contra2/01_APOE4_Strat/06_GNOVA/munged/${file_name}_munged" \
        --a1 reference_allele --a2 other_allele --N-col n_samples --snp rs_number \
        --sumstats "$file" \
        --p P
    echo -e "Finished processing file: $file \n" 
done
popd