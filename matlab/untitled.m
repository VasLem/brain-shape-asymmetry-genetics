close all
modality = 'asymmetry';
reduction=1;
clusterDir = ['../results/' modality '/hierarchicalClusteringDemo/STAGE00DATA/'];
template = load([clusterDir modality '_reduction' num2str(reduction) '/levels4/input_info.mat']).preprocTemplate;

f = figure;
showPaintedDoubleFace(f, template, nan)
