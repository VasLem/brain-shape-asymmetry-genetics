%% testing cart2bary
close all;clear all;clear classes;
load('E:\MatlabProjects\MARKSHRIVER\SingleGeneEffects\DATA\SignGenesResults\MorphAssessmentNESTED_ADAMTS2.mat')
load('E:\MatlabProjects\MARKSHRIVER\SingleGeneEffects\DATA\Science\ColorMorphs\Average_Original.mat')

transferPoseLM(MorphAss.Scan,origscan);

Vertices = origscan.Vertices;
Faces = origscan.Faces;
cart = origscan.PoseLM.Vertices;
Findex = origscan.PoseLM.Fi;

[bar] = cart2bary(Vertices,Faces,cart,Findex);
[cart2] = bary2cart(Vertices,Faces,bar,Findex);

