#!/bin/bash
pushd /data/h_vmac/contra2/01_APOE4_Strat/03_GWAMA/results

for file in *cohorts.out; do
    if [ -f "$file" ]; then
        awk '{ print $1" "$2" "$3" "$5" "$6" "$16" "$4" +" }' "$file" | \
        sed '1s/rs_number/MARKERNAME/g; 1s/reference_allele/EA/g; 1s/other_allele/NEA/g; 1s/n_samples/N/g; 1s/eaf/EAF/g; 1s/+/STRAND/g; 1s/beta/BETA/g; 1s/se/SE/g' > \
        "/data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/input/CrossAnc_$file.GWAMA"
        echo "Processing $file done."
    fi
done
popd

#Move to input directory
pushd /data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/input

#Do this following code in the CrossAnc/input directory with all the GWAMA outputs (after GWAMA conversion) to make .in files 
phenotypes=("MEM" "memslopes" "EXF" "exfslopes" "LAN" "lanslopes")

# Loop through each phenotype
for phenotype in "${phenotypes[@]}"; do
    # Use ls and grep to filter filenames for APOE4neg and store in a list file
    ls | grep "${phenotype}.*APOE4neg.*GWAMA" > "CrossAnc_${phenotype}_APOE4neg.in"  
    # Use ls and grep to filter filenames for APOE4pos and store in a list file
    ls | grep "${phenotype}.*APOE4pos.*GWAMA" > "CrossAnc_${phenotype}_APOE4pos.in"  
    # Use ls and grep to filter filenames for ALL and store in a list file
    ls | grep "${phenotype}.*ALL.*GWAMA" > "CrossAnc_${phenotype}_ALL.in"  
done

ls *.in > CrossAnc.driver
popd

echo "Done"