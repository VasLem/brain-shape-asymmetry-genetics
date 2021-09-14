clear
close all
atlas = 'Desikan_Killiany';
lDesikan = load(['../SAMPLE_DATA/atlasses/L' atlas '_Average.mat']).(atlas);
rDesikan = load(['../SAMPLE_DATA/atlasses/R' atlas '_Average.mat']).(atlas);
mask = load('/usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat').index;
lIndices = lDesikan.index(mask);
rIndices = rDesikan.index(mask);
load ../results/data
