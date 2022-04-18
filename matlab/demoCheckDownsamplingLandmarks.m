
clear, close all
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('BrainAsymmetrySignificanceAnalysis'));
brainSurface = load('../SAMPLE_DATA/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/RH/RENDERMATERIAL.mat');
template = brainSurface.RefScan;
reductions = [1, 0.5, 0.1, 0.05, 0.01];
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
% f = figure;
% for c=1:length(templates)    
%     s = size(templates(c).Vertices, 1);
%     for j=1:2
% 
%     axes = subplot(length(templates),2, 2*(c-1) + j);
%     template =  clone(templates(c));
%     template.Material = 'Dull';
%     template.ViewMode = 'solid';
%     template.Visible = true;
%     template.PatchHandle.FaceColor = 'flat';
%     template.RenderAxes = axes;
%     view(axes, (-1)^j * 90,0);
%     axis(axes,'image');
%     axis(axes,'off');
%     colorbar(axes,'off');
%     light = camlight(axes,'headlight');
%     set(light,'Position',get(axes, 'CameraPosition'));
%     title({['Scale:' num2str(reductions(c))],  [num2str(s) '  landmarks']})
%     pos = get(axes, 'Position');
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
%     end
% end
f = figure;
f.Units = 'normalized';
f.Position = [0.2,0.2,0.5,0.5];
 t = tiledlayout(length(reductions),2);
 for c=1:length(templates)   
 t.TileSpacing = 'none';
t.Padding = 'tight';
tm = templates(c);
tm.Vertices(:,1) = -tm.Vertices(:,1);
ax = showPaintedDoubleFace(f,templates(c),nan,nan, [nexttile(t), nexttile(t)]);
daspect(ax(1), [1 1 1]);
daspect(ax(2), [1 1 1]);
 end
 outdir = '../results/asymmetry/demo_asymmetry/';
 if ~isfolder(outdir), mkdir(outdir); end
saveas(f, [outdir '/reduction.svg'])

