#!/bin/bash
pushd /data/h_vmac/contra2/01_APOE4_Strat/05_MAGMA/annotation_files

# All Races
magma --annotate \
    --snp-loc /data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files/AllRaces/AllRaces_combined.bim \
    --gene-loc /data/h_vmac/Programs/magma/NCBI38.gene.loc.symbols \
    --out /data/h_vmac/contra2/01_APOE4_Strat/05_MAGMA/annotation_files/AllRaces
# NHW
magma --annotate \
    --snp-loc /data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files/NHW/NHW_combined.bim \
    --gene-loc /data/h_vmac/Programs/magma/NCBI38.gene.loc.symbols \
    --out /data/h_vmac/contra2/01_APOE4_Strat/05_MAGMA/annotation_files/NHW

popd