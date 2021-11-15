%% Investigating LEFT - RIGHT asymmetry
close all;clear;
restoredefaultpath;

addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('BrainAsymmetrySignificanceAnalysis'));
%% MAIN PARAMETERS
SEED = 42;
APPLY_ON_ATLAS = false;
N_PICKS = 5;
N_SAMPLES_PER_PICK = 50;
% N_SAMPLES_PER_PICK = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
N_REP = 1; %N_REP=1 for AMMI, N_REP>1 for ANOVA
N_ITER = 10000;
THREADS = 12;
REDUCE = 0.1;
% REDUCE = 1;
GPA_REPS = 3;
PERFORM_EXPERIMENTS = false;
%%

rng(SEED); % For reproducible results
DATA_DIR = '../SAMPLE_DATA/';
OUTPUT_DIR = '../results/demo_asymmetry/';

Ns = floor(1/REDUCE);

PREPROC_DIR = [OUTPUT_DIR 'reduce' num2str(Ns) '_gpaReps' num2str(GPA_REPS) '/'];
ANALYS_DIR = [PREPROC_DIR 'N_REP' num2str(N_REP) '/'];
nSamplesPerPickStr = cellstr(num2str(N_SAMPLES_PER_PICK));
ANOVA_DIR = [ANALYS_DIR  'iters' num2str(N_ITER) '_nPicks' num2str(N_PICKS) '_perPickSize'...
    sprintf('%s_',nSamplesPerPickStr{1:end-1},nSamplesPerPickStr{end})  '/'];

if N_REP > 1
    if ~isfolder(ANOVA_DIR), mkdir(ANOVA_DIR); end
else
    if ~isfolder(ANALYS_DIR), mkdir(ANALYS_DIR); end
end
experimentId = '';
if APPLY_ON_ATLAS
    experimentId = [experimentId '_atlasDK'];
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

%% WITH THE UNCORRECTED DATA SET

outPreprocPath = [PREPROC_DIR, 'preproc_data.mat'];
try
    load(outPreprocPath);
    disp(['Loaded preprocessed data from ' outPreprocPath]);
catch
    disp("Preprocessing Data...")
    DATA = cell(1,2);
    
    for r=1:nR
        regphenopath = [phenopath Regions{r} '/'];
        disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
        DATA{r} = load([regphenopath 'STAGE00DATA']);
    end
    
    template = clone(DATA{2}.Region.AvgShape);
    
    LH = DATA{1}.Region.AlignedShapes;
    RH = DATA{2}.Region.AlignedShapes;
    phenoIID = DATA{1}.Region.IID(1:size(LH,3));
    
    %%
    brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
    refTemplate = brainSurface.RefScan;
    
    
    
    %% Preprocess dataset and the test/retest one
    
    sizeD = size(LH, 3);
    testRetestDataset = loadTestRetestDataset();
    tRShape = size(testRetestDataset.LH);
    tRLH = permute(reshape(permute(testRetestDataset.LH, [3,4, 1,2]), [tRShape(3) * tRShape(4), tRShape(1:2)]), [2, 3, 1]);
    tRRH = permute(reshape(permute(testRetestDataset.RH, [3,4, 1,2]), [tRShape(3) * tRShape(4), tRShape(1:2)]), [2, 3, 1]);
    
    [preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] =...
        preprocessSymmetry(refTemplate, cat(3, LH, tRLH), cat(3, RH, tRRH), [phenoIID;cell(40,1)], REDUCE, 1, GPA_REPS);
    
    alignedTRLH = preprocLH(:, :, sizeD + 1: end);
    alignedTRRH = preprocRH(:, :, sizeD + 1: end);
    preprocLH = preprocLH(:, :, 1: sizeD);
    preprocRH = preprocRH(:, :, 1: sizeD);
    preprocPhenoIID = preprocPhenoIID{1:sizeD};
    
    testRetestVariance = testRetestComputeVariance(alignedTRLH, alignedTRRH, preprocTemplate, PREPROC_DIR, 2);
    clear LH RH DATA brainSurface  
    save([PREPROC_DIR, 'preproc_data.mat'], "preprocRH", "preprocLH", "preprocPhenoIID", "preprocTemplate", "template", "preprocLandmarksIndices", "testRetestVariance", '-v7.3');
end
nLandmarks = length(preprocLandmarksIndices);
nSamples = size(preprocLH,3);
templateAdjacency = template.Adjacency;


repPreprocRH = zeros([size(preprocRH), N_REP + 1], 'single');
repPreprocLH = zeros([size(preprocRH), N_REP + 1], 'single');
repPreprocRH(:, :, :, 1) = single(preprocRH);
repPreprocLH(:, :, :, 1) = single(preprocLH);
if N_REP > 1
    mstdLH = sqrt(testRetestVariance.LH);
    mstdRH = sqrt(testRetestVariance.RH);
    for i = 2: N_REP + 1
        repPreprocLH(:, :, :, i) = single(preprocLH + randn(size(preprocLH)) .* mstdLH);
        repPreprocRH(:, :, :, i) = single(preprocRH + randn(size(preprocRH)) .* mstdRH);
    end
else
    repPreprocRH(:,:,:,2) = preprocRH;
    repPreprocLH(:,:,:,2) = preprocLH;
end
clear preprocRH preprocLH
repPreprocShapes = cat(3, repPreprocRH, repPreprocLH);
clear repPreprocRH repPreprocLH 
repPreprocShapes = permute(repPreprocShapes,[2 1 3, 4]);
repPreprocShapes = permute(reshape(repPreprocShapes, 3 * nLandmarks, (2 * nSamples) ,N_REP + 1), [2, 1, 3]) ;
shapes = repPreprocShapes(:,:,1);
repPreprocShapes = repPreprocShapes(:,:,2:end);
mult = double(intmax('int16') - 1) / (max(abs(repPreprocShapes),[],'all'));
RepShapesInt16 = int16(repPreprocShapes .* mult);
if N_REP > 1
    figure; histogram(reshape(RepShapesInt16 - int16(mult * shapes), 1,[])); title({'Integer Landmark Coordinate dislocation','for generated replications'});
    figure; histogram(reshape(repPreprocShapes - shapes, 1,[])); title({'Landmark Coordinate dislocation','for generated replications'});
end
shapes = permute(reshape(shapes', 3, nLandmarks, 2*nSamples), [2,1,3]);
clear repPreprocShapes
X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:);


%% TWO WAY PROCRUSTES ANOVA ON RANDOM SUBSETS OF THE DATA

atlas = loadAtlas('Desikan_Killiany');
if N_REP == 1
    RESULTS_DIR = ANALYS_DIR;
    [out, scores] = computeAmmiModel(shapes); %650 provided for 0.1
    for comp=1:2
        f=figure;
        ax1 = subplot(1,2,1);
        renderBrainSurface(clone(preprocTemplate), scores(:,comp), ax1);
        view(ax1, -90, 0);
        light = camlight(ax1, 'headlight');
        set(light,'Position', get(ax1,'CameraPosition'));
        ax2 = subplot(1,2,2);
        renderBrainSurface(clone(preprocTemplate), scores(:,comp), ax2);
        view(ax2, 90, 0);
        light = camlight(ax2, 'headlight');
        set(light,'Position', get(ax2,'CameraPosition'));
        
        sgtitle(['PC' num2str(comp)]);
        colorbar(axes, 'SouthOutside');
        savefig(f, [RESULTS_DIR 'ammi_pc' num2str(comp) '.fig']);
        saveas(f, [RESULTS_DIR 'ammi_pc' num2str(comp) '.png']);
    end
    %     ,
    N_SAMPLES_PER_PICK = nSamples;
else
    RESULTS_DIR = ANOVA_DIR;
    disp(['Replication-Based Asymmetry Analysis, using ' num2str(N_PICKS) ' random ' num2str(N_SAMPLES_PER_PICK) ' samples selections out of the original dataset.'])
    if APPLY_ON_ATLAS
        [setOut, avgOut, stdOut] = AsymmetryAnalysisOnSubsets(X1,X2, N_ITER, N_SAMPLES_PER_PICK,N_PICKS, mult,1,atlas.index(preprocLandmarksIndices));
    else
        [setOut, avgOut, stdOut] = AsymmetryAnalysisOnSubsets(X1,X2, N_ITER, N_SAMPLES_PER_PICK,N_PICKS, mult,1); %#ok<UNRCH>
    end
    
    out = avgOut;
    
    if PERFORM_EXPERIMENTS
        [retS, retL, retR, retP] = ProcrustesAnova2WayAsymmetryDebuggingExperiments(X1, X2, mult); %#ok<UNRCH>
    end
    
end
%%
% Upsampling
outu = upsampleAnovaStats(out, templateAdjacency, preprocLandmarksIndices);

showstruct = outu;
showPerm=1;
%%
data = ProcrustesAnova2WayAsymmetryOutputProcess(...
    template, showstruct, N_SAMPLES_PER_PICK , showPerm, [RESULTS_DIR 'data' experimentId '.mat'], 1/N_ITER);

f = visualizeBrainAsymmetryData(data, [RESULTS_DIR 'results' experimentId]);

%%
system('git add *');
message = ['AutoUpdate ' datestr(datetime('now'))];
system(['git commit -m "' message '"']);
system('git push origin');


function variance = testRetestComputeVariance(alignedPLH, alignedPRH, template, resultsDir, n_rep)
arguments
    alignedPLH
    alignedPRH
    template
    resultsDir
    n_rep=2
end

shape = size(alignedPLH);
alignedLH = permute(reshape(permute(alignedPLH, [3, 1, 2]), [shape(3) / n_rep, n_rep, shape(1), shape(2)]), [3,4,1,2]);
alignedRH = permute(reshape(permute(alignedPRH, [3, 1, 2]), [shape(3) / n_rep, n_rep, shape(1), shape(2)]), [3,4,1,2]);

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

