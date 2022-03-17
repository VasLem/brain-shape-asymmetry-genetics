clear
close all
addpath(genpath('AIDFUNCTIONS'));
DATA_DIR = '../SAMPLE_DATA/';
DATASET_INDEX = 1;

switch DATASET_INDEX
    case 1
        UKBIOBANK = 'UKBIOBANK';
        DATASET_NAME = 'STAGE00DATA';
        GENO_ID = 'sel19908';
    case 2
        UKBIOBANK = 'MY_UKBIOBANK';
        DATASET_NAME = 'BATCH2_2021_DATA';
        GENO_ID = 'sel16875_rmrel';
end
in = load([DATA_DIR, '/HumanConnectomeProject/SubcorticalMask_HCP.mat']);

phenopath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/PHENOTYPES/'];
genopath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/GENOTYPES/'];

Regions = {'LH' 'RH'};
nR = length(Regions);
DATA = cell(1,2);

for r=1:nR
    regphenopath = [phenopath Regions{r} '/'];
    disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
    DATA{r} = load([regphenopath DATASET_NAME '.mat']);
end
LH = DATA{1}.Region.AlignedShapes;
RH = DATA{2}.Region.AlignedShapes;
phenoIID = DATA{1}.Region.IID(1:size(LH,3));
Region =  DATA{1}.Region;
brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
refTemplate = brainSurface.RefScan;
[preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, 1, 1, 3);
%% Get average shapes
avgLH = mean(LH,3);
avgRH = mean(RH,3);
%%
templateLH = clone(preprocTemplate);
templateRH = clone(preprocTemplate);
close all
avgTemplate = clone(preprocTemplate);
reflRH = avgRH;
reflRH(:,1) = -reflRH(:,1);
avgTemplate.Vertices = (avgLH + reflRH)/2;
origin = mean(avgTemplate.Vertices);
distOriginRH = vecnorm(reflRH - origin,2,2);
distOriginLH = vecnorm(avgLH - origin,2,2);
diffLR = distOriginLH - distOriginRH;
normDiffLR = 2 *( (diffLR - min(diffLR)) / (max(diffLR) - min(diffLR)) - 0.5);
avgTemplate.VertexValue =normDiffLR;
map = customcolormap_preset('red-white-blue');
avgTemplate.ColorMode = 'Indexed';
f = figure;
f.Units = 'normalized';
f.Position = [0.2,0.2,0.5,0.5];
 t = tiledlayout(1,2);
 t.TileSpacing = 'none';
t.Padding = 'tight';
ax = showPaintedDoubleFace(f,avgTemplate,nan,nan, [nexttile(t), nexttile(t)]);
colormap(map)
daspect(ax(1), [1 1 1]);
daspect(ax(2), [1 1 1]);
h = axes(f,'visible','off');
h.Units = 'normalized';
h.Position = [0.1,0,0.8,0.7];
ticks = linspace(-1,1, 10);
cb = colorbar(h, 'South','TickLabels',  round(ticks*10)/10, ...
        'Ticks',linspace(0,1, 10),'FontWeight','bold', 'FontSize',20);
cb(1).Label.String = '\mid L-C \mid_2 - \mid R - C\mid_2';
saveas(f, sprintf('../results/da_visualization.svg'))