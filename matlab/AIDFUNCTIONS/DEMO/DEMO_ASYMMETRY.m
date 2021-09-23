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
applyOnAtlas = true;
% THREADS = 8;atlas
nPicks = 3;
nSamplesPerPick = 100;
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
[preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, reduce, subsample);
nSamples = size(preprocLH,3);



%%SAMPLE_DATA
% clear DATA;
%%
symmetric = (preprocLH + preprocRH)/2;
asymmetric = (preprocLH - preprocRH);
%%
arr = asymmetric;
arr = reshape(arr, [size(arr,1) * size(arr,2) size(arr,3)])';
[coeff, score, latent, tsquared, explained] = pca(arr, "NumComponents", 100);
f = figure;
plot(cumsum(explained));
savefig(f, [RESULTS_DIR, 'asymmetricPCAExplained.fig']);

%%

numLevels = 3;
clustered = hierarchicalClustering(asymmetric,numLevels,true,3,seed);
fig = paintClusters(clustered, preprocTemplate, numLevels);
savefig(fig, [RESULTS_DIR 'asymmetricClustering.fig'])
saveas(fig, [RESULTS_DIR 'asymmetricClustering.png']);
clustered = hierarchicalClustering(symmetric,numLevels,true,3,seed);
fig = paintClusters(clustered, preprocTemplate, numLevels);
savefig(fig, [RESULTS_DIR 'symmetricClustering.fig'])
saveas(fig, [RESULTS_DIR 'symmetricClustering.png']);
%%

%%
% clear LH  RH;
%%

Ns = floor(1/reduce);
nLandmarks = length(preprocLandmarksIndices);

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

%%
% display3DLandmarksArrows(template, AvgShape);
%%

repPreprocRH = zeros([size(preprocRH), nRep + 1], 'single');
repPreprocLH = zeros([size(preprocRH), nRep + 1], 'single');
repPreprocRH(:, :, :, 1) = single(preprocRH);
repPreprocLH(:, :, :, 1) = single(preprocLH);
if nRep > 1
    variance = load('../results/test_retest_information.mat').variance;
    mstdLH = sqrt(variance.LH(preprocLandmarksIndices, :));
    mstdRH = sqrt(variance.RH(preprocLandmarksIndices, :));
    for i = 2: nRep + 1
        repPreprocRH(:, :, :, i) = single(preprocRH + randn(size(preprocRH)) .* mstdRH) ;
        repPreprocLH(:, :, :, i) = single(preprocLH + randn(size(preprocLH)) .* mstdLH);
    end
else
    repPreprocRH(:,:,:,2) = preprocRH;
    repPreprocLH(:,:,:,2) = preprocLH;
end

templateAdjacency = template.Adjacency;
repPreprocShapes = cat(3, repPreprocRH, repPreprocLH);
repPreprocShapes = permute(repPreprocShapes,[2 1 3, 4]);
repPreprocShapes = permute(reshape(repPreprocShapes, 3 * nLandmarks, (2 * nSamples) ,nRep + 1), [2, 1, 3]) ;
shapes = repPreprocShapes(:,:,1);
repPreprocShapes = repPreprocShapes(:,:,2:end);
mult = double(intmax('int16') - 1) / (max(abs(repPreprocShapes),[],'all'));
RepShapesInt16 = int16(repPreprocShapes.*mult);
if nRep > 1
    figure; histogram(reshape(RepShapesInt16 - int16(mult * shapes), 1,[])); title({'Integer Landmark Coordinate dislocation','for generated replications'});
    figure; histogram(reshape(repPreprocShapes - shapes, 1,[])); title({'Landmark Coordinate dislocation','for generated replications'});
end
shapes = permute(reshape(shapes', 3, nLandmarks, 2*nSamples), [2,1,3]);

X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:);
%% TWO WAY PROCRUSTES ANOVA ON RANDOM SUBSETS OF THE DATA

    
if nRep == 1
    out = computeAmmiModel(shapes);
    nSamplesPerPick = nSamples;
else
    disp(['Replication-Based Asymmetry Analysis, using ' num2str(nPicks) ' random ' num2str(nSamplesPerPick) ' samples selections out of the original dataset.'])
    if applyOnAtlas
        atlas = loadAtlas('Desikan_Killiany');
        [setOut, avgOut,stdOut] = AsymmetryAnalysisOnSubsets(X1,X2,nSamplesPerPick,nPicks, nIter,mult,1,atlas.index(preprocLandmarksIndices));
    else
        [setOut, avgOut,stdOut] = AsymmetryAnalysisOnSubsets(X1,X2,nSamplesPerPick,nPicks, nIter,mult,1); %#ok<UNRCH>
    end

    out = avgOut;
    
    if performExperiments
        [retS, retL, retR, retP] = ProcrustesAnova2WayAsymmetryDebuggingExperiments(X1, X2, mult); %#ok<UNRCH>
    end
end
%%
% Upsampling
outu = upsampleAnovaStats(out, templateAdjacency, preprocLandmarksIndices);

showstruct = outu;
showPerm=1;
data = ProcrustesAnova2WayAsymmetryOutputProcess(...
    template, showstruct, nSamplesPerPick , showPerm, [RESULTS_DIR 'data_' experimentName '.mat']);

f = visualizeBrainAsymmetryData(data, [RESULTS_DIR 'results_' experimentName]);

%%
if nRep == 1
    phenoT = table(preprocPhenoIID, preprocPhenoIID, out.Raw.F');
    writetable(phenoT,[RESULTS_DIR 'fluctuatingAMMI.txt'],'WriteVariableNames',false,'Delimiter',' ');
end

%%
system('git add *');
message = ['AutoUpdate ' datestr(datetime('now'))];
system(['git commit -m "' message '"']);
system('git push origin');


