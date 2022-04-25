%ENV TO SET:
%%%%%%%%%%%%
%DATA_ROOT: the directory of the DATA (../SAMPLE_DATA)
%DATASET_INDEX: dataset to use, 1 or 2 (1)
%RESULTS_ROOT: the directory of the results (../results)
%SCRATCH_ROOT: the directory to use for the intermediate files (../results)
%THREADS: Number of threads to use (auto set by local)
%%%%%%%%%%%%
%% Applying Hierarchical Clustering to Brain Symmetry Related information
close all;clear;
if ~isdeployed
    restoredefaultpath;
    addpath(genpath('AIDFUNCTIONS'));
    addpath(genpath('FUNCTIONS'));
    addpath(genpath('.'));
end
NUM_LEVELS = 4;
MAX_NUM_PCS = 500;
DEFAULT_DATASET_INDEX = 1;
REDUCTION_RATE = 1;
DEFAULT_MODALITY = 'asymmetry'; %asymmetry,symmetry
SEED = 42;

rng(SEED); % For reproducible results
if endsWith(cd, "AIDFUNCTIONS/DEMO")
    cd("../..")
end

DATA_DIR = getenv('DATA_ROOT');
if(isempty(DATA_DIR))
    DATA_DIR = '../SAMPLE_DATA/';
end

DATASET_INDEX = getenv('DATASET_INDEX');
if (isempty(DATASET_INDEX))
    DATASET_INDEX = DEFAULT_DATASET_INDEX;
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
RESULTS_ROOT = getenv('RESULTS_ROOT');
if(isempty(RESULTS_ROOT))
    RESULTS_ROOT = '../results/';
end

SCRATCH_ROOT = getenv('SCRATCH_ROOT');
if(isempty(SCRATCH_ROOT))
    SCRATCH_ROOT = '../results/';
end

MODALITY=getenv('MODALITY');
if isempty(MODALITY)
    MODALITY = DEFAULT_MODALITY;
end

SCRATCH_DIR = [SCRATCH_ROOT, MODALITY '/hierarchicalClusteringDemo/' DATASET_NAME '/'];

RESULTS_DIR = [RESULTS_ROOT, MODALITY '/hierarchicalClusteringDemo/' DATASET_NAME '/'];






THREADS = getenv('THREADS');
if(isempty(THREADS))
    THREADS = parcluster('local');
    THREADS = THREADS.NumWorkers;
else
    if ~isnumeric(THREADS)
        THREADS=str2double(THREADS);
    end
end
%% SETTING UP COMPUTATION POWER
try
    parpool(THREADS);
catch
end


SELECTED_VARIANCE_THRESHOLD = 80;


SELECTION_DIR   = [RESULTS_DIR MODALITY '_reduction' num2str(round(1/REDUCTION_RATE)) '/'];
SELECTION_SCRATCH_DIR   = [SCRATCH_DIR MODALITY '_reduction' num2str(round(1/REDUCTION_RATE)) '/'];
SEGMENTATION_DIR = [SELECTION_DIR 'levels' num2str(NUM_LEVELS) '/'];
SEGMENTATION_SCRATCH_DIR = [SELECTION_SCRATCH_DIR 'levels' num2str(NUM_LEVELS) '/'];
SEGMENTATION_SCRATCH_INPUT_OUT = [SEGMENTATION_SCRATCH_DIR, 'input.mat'];
SEGMENTATION_INPUT_OUT = [SEGMENTATION_DIR, 'input_info.mat'];
SEGMENTATION_INPUT_INFO_PROC = ~isfile(SEGMENTATION_INPUT_OUT);
EXPLAINED_COV_OUT = [SELECTION_DIR  'explained_by_cov_dk.mat'];
EXPLAINED_COV_PROC = ~isfile(EXPLAINED_COV_OUT);
RESIDUALS_OUT = [SELECTION_DIR  'residuals.mat'];
RESIDUALS_PROC = ~isfile(RESIDUALS_OUT);
RV_MATRIX_OUT = [SEGMENTATION_DIR 'similarityMatrix.mat'];
RV_MATRIX_PROC = ~isfile(RV_MATRIX_OUT);
SEGMENTATION_OUT = [SEGMENTATION_DIR 'segmentation.mat'];
SEGMENTATION_ORI_OUT = strrep(SEGMENTATION_OUT, 'BATCH2_2021_DATA', 'STAGE00DATA');
SEGMENTATION_PROC = ~isfile(SEGMENTATION_OUT);
FINAL_RET_INFO_PATH =  [SEGMENTATION_DIR, 'phenotype_info_varThres' num2str(SELECTED_VARIANCE_THRESHOLD) '.mat'];
FINAL_RET_PATH = [SEGMENTATION_DIR, 'phenotype_varThres' num2str(SELECTED_VARIANCE_THRESHOLD) '.mat'];
FINAL_RET_PROC = ~isfile(FINAL_RET_PATH);
FINAL_RET_PROC = 1;
FINAL_PARTITIONS_SIZES_PATH = [SEGMENTATION_DIR, 'phenotype_part_sizes_varThres' num2str(SELECTED_VARIANCE_THRESHOLD) '.mat'];
FINAL_PARTITIONS_INFO_ORIG_PATH = strrep(FINAL_PARTITIONS_SIZES_PATH, 'BATCH2_2021_DATA', 'STAGE00DATA');

SEGMENTATION_SCRATCH_INPUT_PROC = isdeployed || ~isfile(SEGMENTATION_SCRATCH_INPUT_OUT);


PERFORM_STEPS_12 =  RV_MATRIX_PROC || RESIDUALS_PROC || EXPLAINED_COV_PROC;

SEGMENTATION_INPUT_PROC = SEGMENTATION_SCRATCH_INPUT_PROC && PERFORM_STEPS_12;
SEGMENTATION_INPUT_PROC = SEGMENTATION_INPUT_PROC || SEGMENTATION_INPUT_INFO_PROC;



if ~isfolder(SEGMENTATION_DIR), mkdir(SEGMENTATION_DIR); end
if ~isfolder(SEGMENTATION_SCRATCH_DIR), mkdir(SEGMENTATION_SCRATCH_DIR); end
diary([SEGMENTATION_DIR  datestr(now) ' log.txt']);
disp(['Location of data: ', DATA_DIR]);
disp(['Using dataset:', num2str(DATASET_INDEX)])
disp(['Location of results: ', RESULTS_DIR]);
disp(['Location of intermediate files: ', SCRATCH_DIR]);

covGenoPath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];
%% GETTING SOME INFO ON THE BRAIN TEMPLATE



if SEGMENTATION_INPUT_PROC
    in = load([DATA_DIR, 'HumanConnectomeProject/SubcorticalMask_HCP.mat']);

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
    brainSurface = load([phenopath 'LH/RENDERMATERIAL.mat']);
    refTemplate = brainSurface.RefScan;

    % GPA
    [preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, REDUCTION_RATE, 1, 3);


    switch MODALITY
        case 'asymmetry'
            A = (preprocLH - preprocRH);
        case 'symmetry'
            A =  (preprocLH + preprocRH) / 2;
    end
    centroidSizesLH = DATA{1}.Region.CentroidSizes;
    centroidSizesRH = DATA{2}.Region.CentroidSizes;
    regionName = Region.Name;
    clear DATA preprocLH preprocRH LH RH Region
    save(SEGMENTATION_INPUT_OUT, "preprocTemplate","preprocPhenoIID", '-v7.3');
    if ~isdeployed
        save(SEGMENTATION_SCRATCH_INPUT_OUT, 'A',  'centroidSizesLH', 'centroidSizesRH', 'regionName', 'preprocLandmarksIndices', '-v7.3');
    end
end
%%
if ~exist('preprocTemplate', 'var')
    load(SEGMENTATION_INPUT_OUT);
end




%% STEP 1:
% Align covariates with phenotype
load(covGenoPath, "COV");
covariates = COV;
covAssignmentMatrix = (str2double(preprocPhenoIID) == str2double(covariates.IID)');
[covPhenoIndex, covIndex] = find(covAssignmentMatrix);
clear covAssignmentMatrix
preprocPhenoIID = preprocPhenoIID(covPhenoIndex);
if  PERFORM_STEPS_12
    if ~exist('A','var')
        load(SEGMENTATION_SCRATCH_INPUT_OUT)
    end
    [nVertices,~,nSubj] = size(A);
    sym2DMatrix = permute(A,[2 1 3]);
    sym2DMatrix = reshape(sym2DMatrix,nVertices*3,nSubj)';
    clear A
    covData = covariates.DATA(covIndex, :);
    sym2DMatrix = sym2DMatrix(covPhenoIndex, :);
    centroidSizesLH = centroidSizesLH(covPhenoIndex);
    centroidSizesRH = centroidSizesRH(covPhenoIndex);
    avgT = mean(sym2DMatrix,1);
    % Removing components that can be controlled from covariates, keeping only
    % the residuals from the corresponding PLS model
end

if EXPLAINED_COV_PROC
    %% Apply covariates control analysis on DK partitions
    % see how much of the variance of each partition set of landmarks is explained by covariates
    atlasName = 'Desikan_Killiany';
    atlas = loadAtlas(atlasName,'L');
    %%
    atlasIndices = atlas.index(preprocLandmarksIndices);
    atlas3DIndices = repmat(atlasIndices,1,3)';
    atlas3DIndices = atlas3DIndices(:);
    %%
    labels = unique(atlas3DIndices);
    explained3DCov = nan * zeros(length(atlas3DIndices),1);
    for k=1:length(labels)
        i=labels(k);
        atlasMask =  atlas3DIndices == i;
        [~, explained3DCov(atlasMask)] = controlForCovariates([covData, centroidSizesLH(:), centroidSizesRH(:)], sym2DMatrix(:, atlasMask));
    end
    explainedCov = explained3DCov(1:3:length(explained3DCov));
    save(EXPLAINED_COV_OUT, 'explainedCov', '-v7.3');
end
%%

if ~exist('explainedCov', 'var')
    load(EXPLAINED_COV_OUT)
end
f = figure;
view(gca,90,0);

visualizeCovExp = clone(preprocTemplate);
visualizeCovExp.VertexValue = explainedCov;
visualizeCovExp.ColorMode = 'indexed';
visualizeCovExp.Material = 'Dull';
visualizeCovExp.ViewMode = 'solid';

showPaintedDoubleFace(f,visualizeCovExp,nan,nan,nan,nan,[0 1],nan,nan)
colorbar(gca,'SouthOutside');
saveas(f, [SELECTION_DIR, 'explainedDKCovariates.svg']);

%%
if RESIDUALS_PROC
    disp("Fitting PLS model to covariates and removing their effect from the phenotype..");
    [resT, explainedCovWhole] = controlForCovariates([covData,  centroidSizesLH(:), centroidSizesRH(:)], sym2DMatrix);
    resT = repmat(avgT,size(resT,1),1)+resT;
    disp("Saving computed residuals..");
    save(RESIDUALS_OUT,'resT', 'explainedCovWhole', '-v7.3');
end

if DATASET_INDEX == 1
    %% STEP 2: BUILDING RV MATRIX
    if RV_MATRIX_PROC
        disp("Computing RV matrix..")
        if ~exist('resT', 'var')
            load(RESIDUALS_OUT);
            disp(['Loaded computed residuals from ' RESIDUALS_OUT]);
        end

        disp("Using covariates controlled phenotype for segmentation..")
        similarityMat = buildRVmatrixDim(resT, 'cov', 3);

        disp(['Saving RV Matrix to ' SEGMENTATION_DIR 'similarityMatrix.mat ..'])
        save(RV_MATRIX_OUT,'similarityMat','-v7.3');
    end

    %% STEP 3: RUNNING CLUSTERING
    if SEGMENTATION_PROC
        if ~exist('similarityMat', 'var')
            load(RV_MATRIX_OUT);
            disp(['Loaded RV matrix from ' RV_MATRIX_OUT])
        end
        clustered = hierarchicalClustering(similarityMat,NUM_LEVELS,true,3,SEED);
        clusterArray = getClustersAsArray(clustered, NUM_LEVELS);
        save(SEGMENTATION_OUT,'clusterArray','clustered','-v7.3');
        clear similarityMat
    else
        load(SEGMENTATION_OUT);
    end
    %%
    [fig, fig2, handles] = paintClusters(clustered, preprocTemplate, NUM_LEVELS, true, [],nan,colorcube(2));
    saveas(fig, [SEGMENTATION_DIR 'segmentation_circular.png']);
    print(fig2,'-dsvg','-r300',[SEGMENTATION_DIR 'segmentation.svg']);
else
    disp("Loading STAGE00DATA partitions..")
    load(SEGMENTATION_ORI_OUT);
end
%% Compute the PCA Components
% Investigate the number of PCA components to keep for each segment based on the explained Variance
% Get the partition of A that corresponds to the specific cluster with C
% points and then transform it into a N*3 x C array
nSubj = length(preprocPhenoIID);
nVertices = size(clusterArray, 2);
partitions_num = 2^(NUM_LEVELS + 1) - 1;
explainedThresholds = 50:10:90;
clusterExpComponentsNum = zeros(partitions_num, length(explainedThresholds));
means = zeros(partitions_num, 3, nSubj);
templateCenters = zeros(partitions_num,3);
if ~exist('resT', 'var')
    load(RESIDUALS_OUT)
end
reshapedResT = reshape(resT', 3, nVertices, nSubj);
clear resT
disp('Constructing histogram showing the number of components/variance explained for each partition..')
ppb = ParforProgressbar(partitions_num);
explained = cell(partitions_num,1);
pcaScores = cell(partitions_num, 1);
parfor c=1:partitions_num
    try
        k = length(de2bi(c));
        clusterVec = clusterArray(k, :);
        i = c - (2^(k-1) - 1); % cluster index
        selectionMask = clusterVec==i;
        if ~any(selectionMask), continue; end
        ftResT = reshape(reshapedResT(:, selectionMask, :), 3 * sum(selectionMask), nSubj)';
        [~, pcaScores{c}, ~, ~, explained{c}] = pca(ftResT, 'Centered', true);
        explainedAccum = cumsum(explained{c});
        clusterExpComponentsNum(c,:) = sum(explainedAccum < explainedThresholds, 1);
        means(c, :, :) = mean(reshapedResT(:, selectionMask, :), 2);
        templateCenters(c, :) = mean(preprocTemplate.Vertices(selectionMask, :), 1);
        ppb.increment();
    catch ME
        disp(c)
        rethrow(ME)
    end
end
delete(ppb);
fig = figure;
bar(explainedThresholds, clusterExpComponentsNum');
xlabel('Explained Variance Threshold');
ylabel('Number of Components Required');
title('Explained Variance~Number Of Components per Cluster');
savefig(fig, [SEGMENTATION_DIR 'explainedVariance.fig']);
saveas(fig, [SEGMENTATION_DIR 'explainedVariance.png']);
%%
if FINAL_RET_PROC
    if DATASET_INDEX == 1
        disp(['Computing PCA that explains ' num2str(SELECTED_VARIANCE_THRESHOLD) '% of variance for each partition..'])
        
        partitionsSizes = zeros(partitions_num, 1);
        explainedVariance = zeros(partitions_num,1) + SELECTED_VARIANCE_THRESHOLD;
        parfor k = 1:partitions_num
            explainedAccum = cumsum(explained{k});
            partitionsSizes(k) =  find(explainedAccum >= SELECTED_VARIANCE_THRESHOLD, 1 ,'first');
            if MAX_NUM_PCS < partitionsSizes(k)
                explainedVariance(k) = explainedAccum(MAX_NUM_PCS);
                partitionsSizes(k) = MAX_NUM_PCS;
            end
        end
        save(FINAL_PARTITIONS_SIZES_PATH, 'partitionsSizes', 'explainedVariance', '-v7');
    else
        disp(['Loading PCA number that explains ' num2str(SELECTED_VARIANCE_THRESHOLD) '% of variance for each partition of STAGE00DATA..'])

        load(FINAL_PARTITIONS_INFO_ORIG_PATH);
    end

    maxNumComponents = max(clusterExpComponentsNum(: , explainedThresholds == SELECTED_VARIANCE_THRESHOLD));
    ppb = ParforProgressbar(partitions_num);
    parfor k = 1:partitions_num
        scores = pcaScores{k};
        clusterPCAPhenoFeatures{k} = scores(:, 1:partitionsSizes(k));
        ppb.increment();
    end
    delete(ppb)
    PHENO=clusterPCAPhenoFeatures;
    PHENO_IID = preprocPhenoIID;
    save(FINAL_RET_INFO_PATH,'explained', 'templateCenters','means','SELECTED_VARIANCE_THRESHOLD', '-v7.3');
    save(FINAL_RET_PATH, 'PHENO', 'PHENO_IID', '-v7.3');
    
end


%%
% Show the behavior of the means of the clustered regions in 3D
if ~exist('templateCenters', 'var')
    load(FINAL_RET_PATH)
end
localBehaviors = zeros(partitions_num,3, 3);
for i=1:partitions_num
    localBehaviors(i,:, :) = pca(squeeze(means(i, :, :))', 'NumComponents',3);
end

maxlims = max(preprocTemplate.Vertices);
minlims = min(preprocTemplate.Vertices);

fig = figure;
for i=1:2
    ax=subplot(1,2,i);
    ax.Visible = false;
    hold on

    if i == 1
        xlim_min = 1.3 * (maxlims(1) + minlims(1)) / 3;
    else
        xlim_min =minlims(1);
    end
    flag = templateCenters(:,1) >= xlim_min;
    quiver3(ax, templateCenters(flag,1), templateCenters(flag,2), templateCenters(flag,3), localBehaviors(flag, 1, 1), localBehaviors(flag, 2, 1), localBehaviors(flag,3,1), 'LineWidth', 3);
    quiver3(ax, templateCenters(flag,1), templateCenters(flag,2), templateCenters(flag,3), localBehaviors(flag, 1, 2), localBehaviors(flag, 2, 2), localBehaviors(flag,3,2), 'LineWidth', 3);
    quiver3(ax, templateCenters(flag,1), templateCenters(flag,2), templateCenters(flag,3), localBehaviors(flag, 1, 3), localBehaviors(flag, 2, 3), localBehaviors(flag,3,3), 'LineWidth', 3);
    legend('PC1', 'PC2', 'PC3');
    hold off
    brainShape = clone(preprocTemplate);
    brainShape.RenderAxes = gca;
    brainShape.Visible=true;
    % brainShape.ViewMode = 'SolidWireframe';
    brainShape.Material = 'Dull';
    if i == 1
        view(90, 0)
    else
        view(-90,0)
    end
    light = camlight(ax,'headlight');
    set(light,'Position',get(ax,'CameraPosition'));
    xlim(ax, [xlim_min, maxlims(1)]);
    daspect(ax(1), [1 1 1]);
end
saveas(fig, [SEGMENTATION_DIR 'clusterPCA.png']);

diary off;

if ~isdeployed
    close all
end

function clusterArray = getClustersAsArray(clustered, depth)
clusterArray = ones(depth, length(clustered.indices));
for d=0: depth
    clusterArray(d + 1, :) = getClustersAsVec(clustered, d);
end
end

function clusterVec = getClustersAsVec(clustered, level)
clusterVec = ones(1, length(clustered.indices));
if ~isempty(clustered.parts) && level~=0
    part1 = getClustersAsVec(clustered.parts{1}, level - 1);
    part2 = getClustersAsVec(clustered.parts{2}, level - 1) + max(part1);
    clusterVec(clustered.parts{1}.indices) = part1;
    clusterVec(clustered.parts{2}.indices) = part2;
end
end
