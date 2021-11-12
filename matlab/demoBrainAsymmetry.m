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
subsample = 0.1;
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

template = clone(DATA{2}.Region.AvgShape);

LH = DATA{1}.Region.AlignedShapes;
RH = DATA{2}.Region.AlignedShapes;
phenoIID = DATA{1}.Region.IID(1:size(LH,3));

%% 
brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
refTemplate = brainSurface.RefScan;



%%
% We are using only one iteration of Procrustes analysis, so that we can
% compare different datasets, otherwise it would not be possible,
% as we would then have different point of reference
% TODO fix this

sizeD = size(LH, 3);
testRetestDataset = loadTestRetestDataset();
tRShape = size(testRetestDataset.LH);
tRLH = permute(reshape(permute(testRetestDataset.LH, [3,4, 1,2]), [tRShape(3) * tRShape(4), tRShape(1:2)]), [2, 3, 1]);
tRRH = permute(reshape(permute(testRetestDataset.RH, [3,4, 1,2]), [tRShape(3) * tRShape(4), tRShape(1:2)]), [2, 3, 1]);

[preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] =...
    preprocessSymmetry(refTemplate, cat(3, LH, tRLH), cat(3, RH, tRRH), [phenoIID;cell(40,1)], reduce, 1, 3);

alignedTRLH = preprocLH(:, :, sizeD + 1: end);
alignedTRRH = preprocRH(:, :, sizeD + 1: end);
preprocLH = preprocLH(:, :, 1: sizeD);
preprocRH = preprocRH(:, :, 1: sizeD);
preprocPhenoIID = preprocPhenoIID{1:sizeD};
%%
testRetestVariance = testRetestComputeVariance(alignedTRLH, alignedTRRH, preprocTemplate, RESULTS_DIR, 2);


nLandmarks = length(preprocLandmarksIndices);


%%SAMPLE_DATA
% clear DATA;
%%
repPreprocRH = zeros([size(preprocRH), nRep + 1], 'single');
repPreprocLH = zeros([size(preprocRH), nRep + 1], 'single');
repPreprocRH(:, :, :, 1) = single(preprocRH);
repPreprocLH(:, :, :, 1) = single(preprocLH);
if nRep > 1
    mstdLH = sqrt(testRetestVariance.LH);
    mstdRH = sqrt(testRetestVariance.RH);
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


function variance = testRetestComputeVariance(alignedPLH, alignedPRH, template, resultsDir, nRep)
arguments
    alignedPLH
    alignedPRH
    template
    resultsDir
    nRep=2
end
    
shape = size(alignedPLH);
alignedLH = permute(reshape(permute(alignedPLH, [3, 1, 2]), [shape(3) / nRep, nRep, shape(1), shape(2)]), [3,4,1,2]);
alignedRH = permute(reshape(permute(alignedPRH, [3, 1, 2]), [shape(3) / nRep, nRep, shape(1), shape(2)]), [3,4,1,2]);

outLH = mean(var(alignedLH,0,4),3);
outRH = mean(var(alignedRH,0,4),3);
%The following is the reason why we decided to keep both variances and not
%compute the mean out of them. For the test retest dataset, technical
%measurement errors do not seem to be characterized by symmetry on the
%midsaggital plane
f = figure;
hist(outLH - outRH);
title({'Difference between measurement variances of contralateral and symmetrical,', 'as of the midsaggital plane, landmarks'})
saveas(f, [resultsDir 'testRetestVarDiff.png'])
variance.LH = outLH;
variance.RH = outRH;
%%
template = clone(template);
f = figure;
axes = subplot(2, 2, 1);
renderBrainSurface(template, mean(variance.LH, 2), axes);
view(axes,-90,0);
light = camlight(axes,'headlight');
set(light,'Position',get(axes,'CameraPosition'));
title('Left');
colorbar(axes, 'SouthOutside');


axes = subplot(2, 2, 2);
renderBrainSurface(template, mean(variance.LH, 2), axes);
view(axes,90,0);
light = camlight(axes,'headlight');
set(light,'Position',get(axes,'CameraPosition'));
title('Left');
colorbar(axes,'SouthOutside');


axes = subplot(2,2,[3,4]);
title(axes,'Right');
axes = subplot(2, 2, 3);
renderBrainSurface(template, mean(variance.RH, 2), axes);
view(axes,-90,0);
light = camlight(axes,'headlight');
set(light,'Position',get(axes,'CameraPosition'));
title('Right');
colorbar(axes,'SouthOutside');

axes = subplot(2, 2, 4);
renderBrainSurface(template, mean(variance.RH, 2), axes);
view(axes,90,0);
light = camlight(axes,'headlight');
set(light,'Position',get(axes,'CameraPosition'));
title('Right');
colorbar(axes,'SouthOutside');

saveas(f, [resultsDir 'test_retest_variance.png']);
save([resultsDir 'test_retest_information.mat'], "variance","-v7");
end

