% Testing batchAssess
%%
%%
close all;clear all;clear classes;
B = batchAssess;
B.InputFolder = 'C:\MATLAB\Work\Projects\WAReferenceRange\DATA\AssymAnalysis\03MappedAC';
B.NormFolder = 'C:\MATLAB\Work\Projects\WAReferenceRange\DATA\AssymAnalysis\04MirroredAC';
B.OutputFolder = 'C:\MATLAB\Work\Projects\WAReferenceRange\DATA\AssymAnalysis\05AssymAssessAC';
B.Significance = 2;
load('C:\MATLAB\Work\Projects\WAReferenceRange\DATA\AssymAnalysis\MaskIndex.mat')
B.MaskIndex = MaskIndex;
%%
load('C:\MATLAB\Work\Projects\WAReferenceRange\DATA\AMTemplate\AM.mat');
viewer(RefScan)
%%
MaskIndex = RefScan.UserData;
%%
process(B);
%%
close all;clear all;clear classes;
B = batchAssessmentAverager;
B.InputFolder = 'C:\MATLAB\Work\Projects\WAReferenceRange\DATA\AssymAnalysis\05AssymAssessAC';
%%
process(B);
%%
AverageAssessment = B.Average;