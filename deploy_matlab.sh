#!/bin/bash

module purge
module load matlab
mkdir -p deployed
cd matlab
echo 1
mcc -mv demoBrainAsymmetryGenome.m -d ../deployed/
echo 2
mcc -mv demoBrainAsymmetryHierarchicalClustering.m -d ../deployed/
cd ..
