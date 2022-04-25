%% how to render on brain/face surfaces
addpath(genpath('HierarchicalPlotting/'));

loadpath = 'HierarchicalPlotting/';
%% BRAIN
close all

%1) load in the brain render structure variable, I gave this seperatly
Render = load([loadpath 'Brain/RENDERMATERIAL.mat'])

%2) open a figure
f = figure; f.Position = [4  66   1332    433];f.Color = [1 1 1];

%3) get some a value per segment (e.g. avg thickness
SegmentVal = 1:285;                                              %%%%% CHANGE %%%%%
cmap = colormap(parula(256));
range = [min(SegmentVal) max(SegmentVal)];% colorbar range      %%%%% CHANGE %%%%%
crit = mean(SegmentVal);% square if SegmentVal > crit threshold; circle if segmentval < crit
showcolorbar = true;

%4) render
fp = plotInBRAINFigurePanel(f,Render,SegmentVal,range,cmap,crit,showcolorbar);
