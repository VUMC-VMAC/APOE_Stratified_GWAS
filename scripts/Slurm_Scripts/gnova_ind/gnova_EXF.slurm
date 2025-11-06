#!/bin/bash
#SBATCH --mail-user=alex.contreras@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --array=1-65
#SBATCH --job-name=GNOVA_EXF_NHW
#SBATCH --output=/data/h_vmac/contra2/01_APOE4_Strat/logs/06_GNOVA/GNOVA_EXF_%A_%a.out

trait_list=$(cat /data/h_vmac/contra2/01_APOE4_Strat/06_GNOVA/traits/trait_list.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}')

#Activate GNOVA Environment
source /data/h_vmac/contra2/01_APOE4_Strat/06_GNOVA/GNOVA_env/bin/activate

#Run GNOVA
sh /data/h_vmac/contra2/01_APOE4_Strat/scripts/Bash_Scripts/GNOVA/GNOVA_NHW_EXF.sh $trait_list

deactivate