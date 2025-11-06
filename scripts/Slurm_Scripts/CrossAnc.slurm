#!/bin/bash

##The below arguments are necessary for a job array in slurm.

#SBATCH --mail-user=alex.contreras@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=24GB
#SBATCH --array=1-18%18
#SBATCH --job-name=Cross_Ancestry
#SBATCH --output=/data/h_vmac/contra2/01_APOE4_Strat/logs/04_CrossAnc/CrossAnc_%A_%a_final.out

#Load R
load_R

#This has the names of the .in files
infile=$(cat /data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/input/CrossAnc.driver | awk -v line="$SLURM_ARRAY_TASK_ID" '{if (NR == line) print $1}') 

#For using the correct bim file:
datarace="AllRaces"

#Directory where .in files are at and .GWAMA files are at:
infile_dir=/data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/input
outfile_dir=/data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/results
genetic_dir=/data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files

scripts_dir=/data/h_vmac/contra2/01_APOE4_Strat/scripts

echo -e "\n$infile\n##Getting all necessary info from file\n"

pushd $infile_dir

cohortCount=$(wc -l < "$infile")
cat "$infile_dir/$infile"

echo -e "The number of cohorts in this fileset is: $cohortCount \n Now running the GWAMA program on this fileset:"

GWAMA --filelist "$infile_dir/$infile" --output "${outfile_dir}/MetaCrossAnc_${infile}" --quantitative --name_marker MARKERNAME

#SNPs must appear in both ancestries
if [ "$cohortCount" -eq 2 ]; then
    echo "Cohort count is 2, so SNP must appear in at least 2 cohorts"
    awk '{ if(($15>=2) || (NR==1)) print}' ${outfile_dir}/MetaCrossAnc_${infile}.out >  ${outfile_dir}/MetaCrossAnc_${infile}_2of2Ancestry.out
    Rscript ${scripts_dir}/R_Scripts/CrossAnc_man_qq_plots.R ${outfile_dir}/MetaCrossAnc_${infile}_2of2Ancestry.out ${genetic_dir}/${datarace}/${datarace}_combined.bim   
    rm ${outfile_dir}/MetaCrossAnc_${infile}.out
else
    echo "Cohort count is out of the range this script expects."
fi

popd