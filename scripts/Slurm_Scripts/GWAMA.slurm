#!/bin/bash

##The below arguments are necessary for a job array in slurm.

#SBATCH --mail-user=alex.contreras@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --array=1-36%36
#SBATCH --job-name=GWAMA
#SBATCH --output=/data/h_vmac/contra2/01_APOE4_Strat/logs/03_GWAMA/GWAMA_APOE_STRAT_%A_%a_final.out

#Load R with my custom load_R function
load_R

#This has the names of the .in files
infile=$(cat /data/h_vmac/contra2/01_APOE4_Strat/02_GWAS/GWAMA_conversion/fileset_GWAMA_driver.txt | awk -v line="$SLURM_ARRAY_TASK_ID" '{if (NR == line) print $1}') 

#Setting up main directories.
dir=/data/h_vmac/contra2/01_APOE4_Strat
plots_output_dir=/data/h_vmac/contra2/01_APOE4_Strat/03_GWAMA/plots
outfile_dir=/data/h_vmac/contra2/01_APOE4_Strat/03_GWAMA/results
infile_dir=/data/h_vmac/contra2/01_APOE4_Strat/02_GWAS/GWAMA_conversion
scripts_dir=/data/h_vmac/contra2/01_APOE4_Strat/scripts
genetic_dir=/data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files

#For using the correct bim file:
race=`echo $infile | cut -d"_" -f1` #NHW or NHB
if [ "$race" == "NHW" ]; then
  datarace="NHW"
else
  datarace="AllRaces"
fi

pushd $infile_dir
cohortCount=$(wc -l < "$infile")
cat "$infile_dir/$infile"
echo -e "The number of cohorts in this fileset is: $cohortCount \n Now running the GWAMA program on this fileset:"

GWAMA --filelist "$infile_dir/$infile" --output "${outfile_dir}/Meta_${infile}_GWAMAoutput" --quantitative --name_marker MARKERNAME

#Incorporating Derek's point of reference rules
if [ "$cohortCount" -eq 3 ]; then
    echo "Cohort count is 3, so SNP must appear in at least 2 cohorts"
    awk '{ if(($15>=2) || (NR==1)) print}' ${outfile_dir}/Meta_${infile}_GWAMAoutput.out >  ${outfile_dir}/Meta_${infile}_GWAMAoutput_2of3cohorts.out
    Rscript ${scripts_dir}/R_Scripts/GWAMA_man_qq_plots.R ${outfile_dir}/Meta_${infile}_GWAMAoutput_2of3cohorts.out ${genetic_dir}/${datarace}/${datarace}_combined.bim
elif [ "$cohortCount" -eq 8 ]; then
    echo "Cohort count is 8, so SNP must appear in at least 6 cohorts"
    awk '{ if(($15>=6) || (NR==1)) print}' ${outfile_dir}/Meta_${infile}_GWAMAoutput.out >  ${outfile_dir}/Meta_${infile}_GWAMAoutput_6of8cohorts.out
    Rscript ${scripts_dir}/R_Scripts/GWAMA_man_qq_plots.R ${outfile_dir}/Meta_${infile}_GWAMAoutput_6of8cohorts.out ${genetic_dir}/${datarace}/${datarace}_combined.bim
elif [ "$cohortCount" -eq 9 ]; then
    echo "Cohort count is 9, so SNP must appear in at least 6 cohorts" 
    awk '{ if(($15>=6) || (NR==1)) print}' ${outfile_dir}/Meta_${infile}_GWAMAoutput.out >  ${outfile_dir}/Meta_${infile}_GWAMAoutput_6of9cohorts.out
    Rscript ${scripts_dir}/R_Scripts/GWAMA_man_qq_plots.R ${outfile_dir}/Meta_${infile}_GWAMAoutput_6of9cohorts.out ${genetic_dir}/${datarace}/${datarace}_combined.bim       
else
    echo "Cohort count is out of the range this script expects."
fi

popd