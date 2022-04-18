
clear, close all
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('BrainAsymmetrySignificanceAnalysis'));
brainSurface = load('../SAMPLE_DATA/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/RH/RENDERMATERIAL.mat');
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
for c=1:length(templates)
    s = size(templates(c).Vertices, 1);
    for j=1:2
        if j==1
            axes = subplot(length(templates)/3, 3, c);
        end
            template =  clone(templates(c));
            if j==1
                template.Vertices(:,2) = template.Vertices(:,2) - abs(max(template.Vertices(:,2))) - 0.05;
                template.Vertices(:,1) = -template.Vertices(:,1);
            else
                template.Vertices(:,2) = template.Vertices(:,2) + abs(min(template.Vertices(:,2))) + 0.05;
            end
        
            template.Material = 'Dull';
            template.ViewMode = 'solid';
            template.Visible = true;
            template.PatchHandle.FaceColor = 'flat';
            template.RenderAxes = axes;
            
            view(axes,  90,0);
            if j==1
            axis(axes,'image');
            axis(axes,'off');
            
            colorbar(axes,'off');
            light = camlight(axes,'headlight');
            set(light,'Position',get(axes, 'CameraPosition'));
            end
            title({['Scale:' num2str(reductions(c))],  [num2str(s) '  landmarks']})
            pos = get(axes, 'Position');
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end
    end
    % f = figure;
    % f.Units = 'normalized';
    % f.Position = [0.1,0.1, 0.8*(1/length(templates)), 0.8];
    %  t = tiledlayout(length(reductions),2);
    %  for c=1:length(templates)
    %  t.TileSpacing = 'none';
    % t.Padding = 'tight';
    % tm = templates(c);
    % tm.Vertices(:,1) = -tm.Vertices(:,1);
    % ax = showPaintedDoubleFace(f,templates(c),nan,nan, [nexttile(t), nexttile(t)]);
    % daspect(ax(1), [1 1 1]);
    % daspect(ax(2), [1 1 1]);
    %  end
    outdir = '../results/downsampling/';
    if ~isfolder(outdir), mkdir(outdir); end
    saveas(f, [outdir 'reduction.svg'])

