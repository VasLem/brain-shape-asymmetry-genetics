%% Investigating LEFT - RIGHT asymmetry
close all;clear;
restoredefaultpath;
%%
addpath(genpath('AIDFUNCTIONS'));
%%
loadWhileInLab = 1;

seed = 42;
rng(seed); % For reproducible results
if endsWith(cd, "AIDFUNCTIONS/DEMO")
    cd("../..")
end
DATA_DIR = '../SAMPLE_DATA/';
RESULTS_DIR = '../results/demo_asymmetry/';
applyOnAtlas = false;
% THREADS = 8;atlas
nPicks = 5;
nSamplesPerPick = 50;
% nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
% reduce = 0.05;
nRep = 3;
nIter = 1000;
THREADS = 8;
reduce = 0.1;
subsample = 1;
% Define following when nRep>1 -> no use of AMMI
% nIter = 5000;
%nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
%nPicks = 10;
performExperiments = 0;
% addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));

Ns = floor(1/reduce);
if loadWhileInLab
    nSamples = 19644;
else
    nSamples = 1000; %#ok<UNRCH>
end

if nRep > 1
    experimentName = [num2str(nSamples) '_' num2str(Ns) '_' ...
        num2str(nIter)  '_' num2str(nPicks) '_' num2str(length(nSamplesPerPick)) '_' num2str(nSamplesPerPick(1))...
        '_' num2str(nSamplesPerPick (end)) '_' num2str(nRep)];
else
    experimentName = [num2str(nSamples) '_' num2str(Ns) ...
        '_' num2str(nRep)];
end
if applyOnAtlas
    experimentName = [experimentName '_DK'];
end
%% SETTING UP COMPUTATION POWER
try
    parpool('local',THREADS);
    %parpool('local',8);
catch
end
%% GETTING SOME INFO ON THE BRAIN TEMPLATE


in = load([DATA_DIR, 'IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat']);

phenopath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/'];
genopath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/'];
covGenoPhenoPath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/COVARIATES/'];

MASK = in.index;
nVertices = length(MASK);
Regions = {'LH' 'RH'};
nR = length(Regions);


%% Extract mat file of serialized objects
% strDATA = cell(1,2);
% strRender = cell(1, 2);
% savePath = [DATA_DIR 'asym_data.mat'];
% disp(['Extracting mat file of serialized objects to ' savePath]);
% for r=1:nR
%     %r=2
%     regphenopath = [phenopath Regions{r} '/'];
%     strDATA{r} =  load([regphenopath 'STAGE00DATA']);
%     strDATA{r}.Region.AvgShape = strDATA{r}.Region.AvgShape.obj2struc();
%     strRender{r} = load([regphenopath 'RENDERMATERIAL.mat']);
%     strRender{r}.RefScan = strRender{r}.RefScan.obj2struc();
%
% end
% save(savePath,  'strDATA', 'strRender');
% clear strDATA strRender;

%% WITH THE UNCORRECTED DATA SET

% phenopath = '/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/';


DATA = cell(1,2);

for r=1:nR
    regphenopath = [phenopath Regions{r} '/'];
    disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
    if ~loadWhileInLab
        DATA{r} = load([regphenopath 'STAGEOODATA']);
        DATA{r}.Region = DATA{r}.Region.Region;
    else
        DATA{r} = load([regphenopath 'STAGE00DATA']);
    end
end

template = clone(DATA{1}.Region.AvgShape);

LH = DATA{1}.Region.AlignedShapes;
RH = DATA{2}.Region.AlignedShapes;
phenoIID = DATA{1}.Region.IID(1:size(LH,3));

%% 
brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
refTemplate = brainSurface.RefScan;



%%
[preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, 1, subsample);

nLandmarks = length(preprocLandmarksIndices);


%%SAMPLE_DATA
% clear DATA;
%%
symmetric = (preprocLH + preprocRH)/2;
asymmetric = (preprocLH - preprocRH);

numLevels = 3;
clustered = hierarchicalClustering(asymmetric,numLevels,true,3,seed);
fig = paintClusters(clustered, preprocTemplate, numLevels);
savefig(fig, [RESULTS_DIR 'asymmetricClustering.fig'])
saveas(fig, [RESULTS_DIR 'asymmetricClustering.png']);
clustered = hierarchicalClustering(symmetric,numLevels,true,3,seed);
fig = paintClusters(clustered, preprocTemplate, numLevels);
savefig(fig, [RESULTS_DIR 'symmetricClustering.fig'])
saveas(fig, [RESULTS_DIR 'symmetricClustering.png']);
