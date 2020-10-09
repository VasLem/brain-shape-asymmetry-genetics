%test kdtree
close all; clear all;clear classes;
load('C:\MATLAB\Work\Projects\WAReferenceRange\AMv12.mat');
load('C:\MATLAB\Work\Projects\WAReferenceRange\DATA\03MappedAC\FSAC_0001_N.mat');
tic;
[Cp,d,I] = kdtree(obj,RefScan);
toc;
tic;
   [N,D] = KNN(RefScan,obj,1);
toc;
