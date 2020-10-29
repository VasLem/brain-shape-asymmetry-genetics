%% Developping a JanusMorpheus model for faces and genes
close all;clear all;clear classes;
StudyPath = 'D:\==Matlab==\MARKSHRIVER\SingleGeneEffects';
cd(StudyPath); cd DATA;
load ABC_ZERO_DFC_LATEST_TO_USE;
load('D:\==Matlab==\MARKSHRIVER\SingleGeneEffects\DATA\NESTED_GENELoopingResult.mat');
RefScan = clone(SymModel.Average);
%%
test = DNA2FaceModel;