#!/bin/bash

module purge
module load matlab/R2019a
mkdir -p deployed
cd deployed

echo 1
mcc -mv demoBrainAsymmetryGenome.m -d ../deployed/
echo 2
mcc -mv demoBrainAsymmetryHierarchicalClustering.m -d ../deployed/
cd ..
