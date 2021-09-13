%% Investigating LEFT - RIGHT asymmetry
close all;clear;
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
% THREADS = 8;
% samplesIndices = 1:1000;
nPicks = 10;
nSamplesPerPick = [200];
% nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
% reduce = 0.05;
% reduce = 0.01;
nRep = 1;
nIter = 1000;
THREADS = 8;
reduce = 0.1;
if loadWhileInLab
    subsample = 0.1;
else
    subsample = 1;
end
% Define following when nRep>1 -> no use of AMMI
% nIter = 5000;
%nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
%nPicks = 10;
performExperiments = 0;
restoredefaultpath;
% addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));

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

Template = clone(DATA{1}.Region.AvgShape);

LH = DATA{1}.Region.AlignedShapes;
RH = DATA{2}.Region.AlignedShapes;
phenoIID = DATA{1}.Region.IID(1:size(LH,3));

%%
Render = cell(1,2);
for r=1:nR
    Render{r} = load([regphenopath 'RENDERMATERIAL.mat']);
end
brainSurface = Render{1};




%%
[preprocTemplate, preprocLH, preprocRH, preprocLandmarksIndices] = preprocessSymmetry(brainSurface.RefScan, LH, RH, reduce, subsample);

%%
assert(min(size(LH,3), sum(str2double(DATA{1}.Region.IID) == str2double(DATA{2}.Region.IID))) == size(LH,3))

%%
checkTemplate =  shape3D;
checkTemplate.Vertices = preprocLH(:,:,1);
checkTemplate.Faces = preprocTemplate.Faces;
checkTemplate.Material = 'Dull';
checkTemplate.ViewMode = 'solid';
checkTemplate.PatchHandle.FaceColor = 'flat';
viewer = checkTemplate.viewer;
light = camlight(viewer.RenderAxes,'headlight');
set(light,'Position',get(viewer.RenderAxes,'CameraPosition'));

%%SAMPLE_DATA
clear DATA;
%%
nSamples = size(preprocLH,3);
%%
symmetric = (preprocLH + preprocRH)/2;
asymmetric = (preprocLH - preprocRH);
numLevels = 3;
clustered = hierarchicalClustering(asymmetric,numLevels,true,3,seed);
fig = paintClusters(clustered, preprocTemplate, numLevels);
savefig(fig, '../results/asymmetricClusterting.fig')
clustered = hierarchicalClustering(symmetric,numLevels,true,3,seed);
fig = paintClusters(clustered, preprocTemplate, numLevels);
savefig(fig, '../results/symmetricClusterting.fig')
%%

%%
clear LH  RH;
%%

Ns = floor(1/reduce);
nLandmarks = length(preprocLandmarksIndices);

if nRep > 1
    experimentName = [num2str(length(samplesIndices)) '_' num2str(Ns) '_' ...
        num2str(nIter)  '_' num2str(nPicks) '_' num2str(length(nSamplesPerPick)) '_' num2str(nSamplesPerPick(1))...
        '_' num2str(nSamplesPerPick (end)) '_' num2str(nRep)];
else
    experimentName = [num2str(length(samplesIndices)) '_' num2str(Ns) ...
        '_' num2str(nRep)];
end

%%
% display3DLandmarksArrows(Template, AvgShape);
%%

repPreprocRH = zeros([size(preprocRH), nRep + 1], 'single');
repPreprocLH = zeros([size(preprocRH), nRep + 1], 'single');
repPreprocRH(:, :, :, 1) = single(preprocRH);
repPreprocLH(:, :, :, 1) = single(preprocLH);

if nRep > 1
    variance = load('test_retest_information.mat').variance;
    rvarLH = variance.LH(preprocLandmarksIndices, :);
    variance.RH(preprocLandmarksIndices, :);
    rvarRH = variance.LH(preprocLandmarksIndices, :);
    for i = 2: nRep + 1
        repPreprocRH(:, :, :, i) = single(preprocRH) + single(randn(size(preprocRH)) .* rvarRH);
        repPreprocLH(:, :, :, i) = single(preprocLH) + single(randn(size(preprocLH)) .* rvarLH);
    end
else
    repPreprocRH(:,:,:,2) = preprocRH;
    repPreprocLH(:,:,:,2) = preprocLH;
end
%%

reducedTemplateAdjacency = Template.Adjacency;
repPreprocShapes = cat(3, repPreprocRH, repPreprocLH);
repPreprocShapes = permute(repPreprocShapes,[2 1 3, 4]);
repPreprocShapes = permute(reshape(repPreprocShapes, 3 * nLandmarks, (2 * nSamples) ,nRep + 1), [2, 1, 3]) ;
shapes = repPreprocShapes(:,:,1);
repPreprocShapes = repPreprocShapes(:,:,2:end);
mult = double(intmax('int16')) / (max(abs(repPreprocShapes),[],'all'));
RepShapesInt16 = int16(repPreprocShapes.*mult);
if nRep > 1
    figure; histogram(reshape(RepShapesInt16 - int16(mult * shapes), 1,[])); title({'Landmark Coordinate dislocation','for generated replications'});
end
%%
X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:);
%%


%% TWO WAY PROCRUSTES ANOVA ON RANDOM SUBSETS OF THE DATA
if nRep == 1
    out = computeAmmiModel(preprocShapes);
    nSamplesPerPick = nSamples;
else
    disp(['Replication-Based Asymmetry Analysis, using ' nPicks ' random ' nSamplesPerPick ' samples selections out of the original dataset.'])
    [setOut, avgOut,stdOut] = AsymmetryAnalysisOnSubsets(X1,X2,nSamplesPerPick,nPicks, nIter,mult,1);
    out = avgOut;
    
    if performExperiments
        [retS, retL, retR, retP] = ProcrustesAnova2WayAsymmetryDebuggingExperiments(X1, X2, mult); %#ok<UNRCH>
    end
end
%%
% Upsampling
outu = upsampleAnovaStats(out, reducedTemplateAdjacency, landmarksIndices);

showstruct = outu;
showPerm=1;
data = ProcrustesAnova2WayAsymmetryOutputProcess(...
    brainSurface, showstruct, nSamplesPerPick , showPerm, [RESULTS_DIR 'data_' experimentName '.mat']);

%%
f = visualizeBrainAsymmetryData(data,[RESULTS_DIR 'results_' experimentName]);

%%
if nRep == 1
    phenoT = table(phenoIID, phenoIID, out.Raw.F');
    writetable(phenoT,[RESULTS_DIR 'fluctuatingAMMI.txt'],'WriteVariableNames',false,'Delimiter',' ');
end

%%
system('git add *');
message = ['AutoUpdate ' datestr(datetime('now'))];
system(['git commit -m "' message '"']);
system('git push origin');


