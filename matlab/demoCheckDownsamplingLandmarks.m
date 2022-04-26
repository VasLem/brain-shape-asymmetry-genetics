
clear, close all
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('BrainAsymmetrySignificanceAnalysis'));
brainSurface = load('../SAMPLE_DATA/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/LH/RENDERMATERIAL.mat');
template = brainSurface.RefScan;
reductions = [1, 0.5, 0.1, 0.05, 0.01, 0.005];
for c = 1: length(reductions)
    if reductions(c) == 1
        templates(c) = template;
        continue;
    end
    disp(['Scale down:', num2str(reductions(c))]);
    [landmarksIndices, reducedFaces, ~]  = getDownsampledLandmarksIndices(template,reductions(c),true);
    reducedTemplate = shape3D;
    reducedTemplate.Vertices = template.Vertices(landmarksIndices, :) ;
    reducedTemplate.Faces = reducedFaces;
    templates(c) = reducedTemplate;
end
%%
f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.1, 0.8/2, 0.8*(2/length(templates))];
p=nan;
for c=1:length(templates)
    s = size(templates(c).Vertices, 1);
    axes = subplot(length(templates)/3, 3, c);
    [ax,p] = showPaintedDoubleFace(f,templates(c),nan,nan, axes,false,[0,1],'white',p);
    title(ax, {['Scale:' num2str(reductions(c))],  [num2str(s) '  landmarks']})
end
outdir = '../results/downsampling/';
set(f, 'InvertHardCopy', 'off');
set(f, 'Color', 'white');
if ~isfolder(outdir), mkdir(outdir); end
saveas(f, [outdir 'reduction.svg'])

