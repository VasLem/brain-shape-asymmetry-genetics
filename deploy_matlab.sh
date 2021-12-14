#!/bin/bash

module purge
module load matlab/R2019a
mkdir -p deployed
cd matlab
echo 1
mcc -mv demoBrainAsymmetryGenome.m -d ../deployed/  -a ./AIDFUNCTIONS -a ./FUNCTIONS -a ./BrainAsymmetrySignificanceAnalysis
echo 2
mcc -mv demoBrainAsymmetryHierarchicalClustering.m -d ../deployed/  -a ./AIDFUNCTIONS -a ./FUNCTIONS -a ./BrainAsymmetrySignificanceAnalysis
cd ..
