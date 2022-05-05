% Create DK atlas image
clear
close all
setuplatex
addpath(genpath('AIDFUNCTIONS'));
atlasL = loadAtlas('Desikan_Killiany','L');
atlasR = loadAtlas('Desikan_Killiany','R');
DATA_DIR = getenv('DATA_ROOT');
if(isempty(DATA_DIR))
    DATA_DIR = '../SAMPLE_DATA/';
end

cmap = brewermap(length(unique(atlasL.index)), 'Spectral');

DATASET_INDEX = getenv('DATASET_INDEX');
if (isempty(DATASET_INDEX))
    DATASET_INDEX = 1;
else
    disp(DATASET_INDEX)
    if ~isnumeric(DATASET_INDEX)
        DATASET_INDEX = str2double(DATASET_INDEX);
    end
end

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
    brainSurfaces{r} = load([regphenopath 'RENDERMATERIAL.mat']);
end

lTemplate = brainSurfaces{1}.RefScan;
rTemplate = brainSurfaces{2}.RefScan;
%%
lTemplate.VertexRGB = zeros(numel(atlasL.index), 3);
exind = atlasL.index(atlasL.index~=-1);
lTemplate.VertexRGB(atlasL.index~=-1, :) = cmap(exind, :);
lTemplate.ColorMode = 'texture';
fig = figure();
fig.Position = [100 100 400 200];
t = tiledlayout(1,1);
t.TileSpacing = 'none';
t.Padding = 'tight';
ax = showPaintedDoubleFace(fig, lTemplate,nan, nan,nexttile(t) , cmap);

saveas(fig, '../results/dk_left.svg')
%%