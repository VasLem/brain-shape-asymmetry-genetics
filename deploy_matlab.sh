#!/bin/bash

module purge
module load matlab/R2019a
mkdir -p deployed
cd matlab
echo 1
${MCC_PATH}mcc -mv -R -singleCompThread -R -nodisplay  demoBrainAsymmetryGenome.m -d ../deployed/ \
 -a ./AIDFUNCTIONS -a ./FUNCTIONS -a ./BrainAsymmetrySignificanceAnalysis \
 -a ../SNPLIB 
echo 2
${MCC_PATH}mcc -mv -R -singleCompThread -R -nodisplay  demoBrainAsymmetryHierarchicalClustering.m -d ../deployed/ \
 -a ./AIDFUNCTIONS -a ./FUNCTIONS -a ./BrainAsymmetrySignificanceAnalysis
cd ..
