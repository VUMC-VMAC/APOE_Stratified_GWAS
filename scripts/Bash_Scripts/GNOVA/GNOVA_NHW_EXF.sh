#!/bin/bash

#Set up array of summary stats
trait_list=$1

cross=EXF
long=exfslopes

input_dir=/data/h_vmac/contra2/01_APOE4_Strat/06_GNOVA/munged
output_dir=/data/h_vmac/contra2/01_APOE4_Strat/06_GNOVA/results
script_dir=/data/h_vmac/contra2/01_APOE4_Strat/scripts/Python_Scripts/GNOVA-2.0

######################################################################################################################################################
## Cross-sectional
######################################################################################################################################################

##Sets up a GNOVA job array for all traits (neg)
python2 $script_dir/gnova.py /data/h_vmac/eissmajm/SUMSTATS/Munged_Traits/${trait_list}.sumstats.gz \
$input_dir/Meta_NHW_${cross}_APOE4neg.in_GWAMAoutput_*_munged.sumstats.gz \
--bfile /data/h_vmac/eissmajm/GNOVA-2.0/genotype_1KG_eur_SNPmaf5/eur_chr@_SNPmaf5 \
--out $output_dir/${cross}/${trait_list}_NHW_${cross}_APOE4neg.txt

##Sets up a GNOVA job array for all traits (pos)
python2 $script_dir/gnova.py /data/h_vmac/eissmajm/SUMSTATS/Munged_Traits/${trait_list}.sumstats.gz \
$input_dir/Meta_NHW_${cross}_APOE4pos.in_GWAMAoutput_*_munged.sumstats.gz \
--bfile /data/h_vmac/eissmajm/GNOVA-2.0/genotype_1KG_eur_SNPmaf5/eur_chr@_SNPmaf5 \
--out $output_dir/${cross}/${trait_list}_NHW_${cross}_APOE4pos.txt

######################################################################################################################################################
## Longitudinal
######################################################################################################################################################
###Sets up a GNOVA job array for all slopes (neg)
python2 $script_dir/gnova.py /data/h_vmac/eissmajm/SUMSTATS/Munged_Traits/${trait_list}.sumstats.gz \
$input_dir/Meta_NHW_${long}_APOE4neg.in_GWAMAoutput_*_munged.sumstats.gz \
--bfile /data/h_vmac/eissmajm/GNOVA-2.0/genotype_1KG_eur_SNPmaf5/eur_chr@_SNPmaf5 \
--out $output_dir/${long}/${trait_list}_NHW_${long}_APOE4neg.txt

#####Sets up a GNOVA job array for all slopes (pos)
python2 $script_dir/gnova.py /data/h_vmac/eissmajm/SUMSTATS/Munged_Traits/${trait_list}.sumstats.gz \
$input_dir/Meta_NHW_${long}_APOE4pos.in_GWAMAoutput_*_munged.sumstats.gz \
--bfile /data/h_vmac/eissmajm/GNOVA-2.0/genotype_1KG_eur_SNPmaf5/eur_chr@_SNPmaf5 \
--out $output_dir/${long}/${trait_list}_NHW_${long}_APOE4pos.txt