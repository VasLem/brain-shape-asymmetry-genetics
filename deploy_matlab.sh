#!/bin/bash
mkdir -p deployed
cd matlab
echo 1
/SOFTWARE/Matlab/R2019a/bin/mcc -mv demoBrainAsymmetryGenome.m -d ../deployed/
echo 2
/SOFTWARE/Matlab/R2019a/bin/mcc -mv demoBrainAsymmetryHierarchicalClustering.m -d ../deployed/
cd ..
