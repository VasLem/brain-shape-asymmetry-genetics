%% Applying Hierarchical Clustering to Brain Symmetry Related information
close all;clear;
restoredefaultpath;
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('FUNCTIONS'));

SEED = 42;
rng(SEED); % For reproducible results
if endsWith(cd, "AIDFUNCTIONS/DEMO")
    cd("../..")
end
DATA_DIR = '../SAMPLE_DATA/';
RESULTS_DIR = '../results/hierarchicalClusteringDemo/';


THREADS = 16;



SELECTION = 'asymmetry';
% SELECTION = 'symmetry';

COV_REMOVAL = 'ccPriorSegmentation';
% COV_REMOVAL = 'ccPostSegmentation';

NUM_LEVELS = 4;
SEGMENTATION_USED = 'mine';
% SEGMENTATION_USED = 'labs';

REDUCTION_RATE = 1;

SELECTION_DIR   = [RESULTS_DIR SELECTION '_reduction' num2str(round(1/REDUCTION_RATE)) '/'];
COV_REMOVAL_DIR = [SELECTION_DIR COV_REMOVAL '/'];
SEGMENTATION_DIR = [COV_REMOVAL_DIR 'levels' num2str(NUM_LEVELS) '_' SEGMENTATION_USED '/'];


if ~isfolder(SEGMENTATION_DIR), mkdir(SEGMENTATION_DIR); end

%% SETTING UP COMPUTATION POWER
try
    parpool('local',THREADS);
catch
end

covGenoPath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/COVARIATES/COVDATAINLIERS.mat'];
%% GETTING SOME INFO ON THE BRAIN TEMPLATE
try
    load([SELECTION_DIR, 'input.mat']);
catch
    
    
    in = load([DATA_DIR, 'IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat']);
    
    phenopath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/'];
    genopath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/'];
    
    Regions = {'LH' 'RH'};
    nR = length(Regions);
    
    DATA = cell(1,2);
    
    for r=1:nR
        regphenopath = [phenopath Regions{r} '/'];
        disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
        DATA{r} = load([regphenopath 'STAGE00DATA']);
    end
    LH = DATA{1}.Region.AlignedShapes;
    RH = DATA{2}.Region.AlignedShapes;
    phenoIID = DATA{1}.Region.IID(1:size(LH,3));
    Region =  DATA{1}.Region;
    brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
    refTemplate = brainSurface.RefScan;
    
    % GPA
    [preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, REDUCTION_RATE, 1, 3);
    
    
    switch SELECTION
        case 'asymmetry'
            A = (preprocLH - preprocRH);
        case 'symmetry'
            A =  (preprocLH + preprocRH)/2;
    end
    clear preprocLH preprocRH LH RH
    centroidSizes = Region.CentroidSizes;
    regionName = Region.Name;
    save([SELECTION_DIR, 'input.mat'], 'A', 'preprocTemplate', 'centroidSizes', 'regionName', 'preprocPhenoIID', 'preprocLandmarksIndices', '-v7.3');
end
%%
preprocTemplate.Vertices(:,1) = -preprocTemplate.Vertices(:,1);
%%

nLandmarks = length(preprocLandmarksIndices);
%% STEP 1:
[nVertices,DIM,nSubj] = size(A);
sym2DMatrix = permute(A,[2 1 3]);
sym2DMatrix = reshape(sym2DMatrix,nVertices*DIM,nSubj)';
clear A
% Align covariates with phenotype
covariates = load(covGenoPath).COV;

%%
covAssignmentMatrix = (str2double(preprocPhenoIID) == str2double(covariates.IID)');
[covPhenoIndex, covIndex] = find(covAssignmentMatrix);

covData = covariates.DATA(covIndex, :);
sym2DMatrix = sym2DMatrix(covPhenoIndex, :);
preprocPhenoIID = preprocPhenoIID(covPhenoIndex);

avgT = mean(sym2DMatrix,1);
% Removing components that can be controlled from covariates, keeping only
% the residuals from the corresponding PLS model
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
    [~, explained3DCov(atlasMask)] = controlForCovariates([covData, centroidSizes(:)], sym2DMatrix(:, atlasMask));
end
%%
explainedCov = explained3DCov(1:3:length(explained3DCov));
f = figure;
view(gca,90,0);
colormap(gca,'jet');
light = camlight(gca,'headlight');
set(light,'Position',get(gca,'CameraPosition'));
drawnow;
visualizeCovExp = clone(preprocTemplate);
visualizeCovExp.VertexValue = explainedCov;

visualizeCovExp.ColorMode = 'indexed';
visualizeCovExp.Material = 'Dull';
visualizeCovExp.ViewMode = 'solid';
colorbar(gca,'SouthOutside');
visualizeCovExp.RenderAxes = gca;
visualizeCovExp.Visible = true;
visualizeCovExp.PatchHandle.FaceColor = 'flat';
axis(gca,'image');
axis(gca,'off');
drawnow;
savefig(f, [SELECTION_DIR, 'explainedDKCovariatesMedial.fig']);
saveas(f, [SELECTION_DIR, 'explainedDKCovariatesMedial.png']);
view(gca,-90,0);
set(light,'Position',get(gca,'CameraPosition'));
savefig(f, [SELECTION_DIR, 'explainedDKCovariatesLateral.fig']);
saveas(f, [SELECTION_DIR, 'explainedDKCovariatesLateral.png']);

%%
diary([SEGMENTATION_DIR  datestr(now) ' log.txt']);
try
    load([SELECTION_DIR  'residuals.mat']);
    disp(['Loaded computed residuals from ' SELECTION_DIR  'residuals.mat']);
catch
    disp("Fitting PLS model to covariates and removing their effect from the phenotype..");
    [resT, explainedCovWhole] = controlForCovariates([covData centroidSizes(:)], sym2DMatrix);
    resT = repmat(avgT,size(resT,1),1)+resT;
    disp("Saving computed residuals..");
    save([SELECTION_DIR  'residuals.mat'],'resT', 'explainedCovWhole', '-v7.3');
end
%% STEP 2: BUILDING RV MATRIX
try
    load([COV_REMOVAL_DIR 'similarityMatrix.mat']);
    disp(['Loaded RV matrix from ' COV_REMOVAL_DIR 'similarityMatrix.mat'])
catch
    disp("Computing RV matrix..")
    if strcmp(COV_REMOVAL, 'ccPriorSegmentation')
        disp("Using covariates controlled phenotype for segmentation..")
        similarityMat = buildRVmatrixDim(resT, 'cov', 3);
    else
        disp("Using covariates uncontrolled phenotype for segmentation..")
        similarityMat = buildRVmatrixDim(sym2DMatrix, 'cov', 3);
    end
    disp(['Saving RV Matrix to ' COV_REMOVAL_DIR 'similarityMatrix.mat ..'])
    save([COV_REMOVAL_DIR 'similarityMatrix.mat'],'similarityMat','-v7.3');
end
%% STEP 3: RUNNING CLUSTERING

if strcmp(SEGMENTATION_USED, 'labs')
    type = 'weiss';%'symmetric laplacian'; %'ratiocute';%'ncute';%'weiss'; is not so 'symmetric laplacian'
    %nice for RH and SH;
    %type = 'weiss';%'ratiocute';%'ncute';%'weiss';
    minPercValue = 1;
    runs = 50;
    [LABELS,MASK] = HierarchicalFacialSegmentationv4(double(similarityMat),NUM_LEVELS,type,runs,minPercValue);
    save([SEGMENTATION_DIR 'segmentation.mat'],'LABELS','MASK','-v7.3');
    
    % VISUALIZING THE SEGMENTATION
    % find the deepest still meaningfull level
    levels = nansum(LABELS,2);
    index = find(levels);
    u_levels = index(end);% levels to analyze
    ULABELS = LABELS(index,:);
    UHI = HierarchicalInterface;
    UHI.nL = u_levels;
    UMASK = MASK(1:UHI.nLC);
    disp(['Number of total clusters: ' num2str(length(find(UMASK)))]);
    v_levels = NUM_LEVELS;% ter visualisatie
    VLABELS = ULABELS(1:v_levels,:);
    VHI = HierarchicalInterface;
    VHI.nL = v_levels;
    VMASK = UMASK(1:VHI.nLC);
    disp(['Number of visualized clusters: ' num2str(length(find(VMASK)))]);
    
    % make figdir
    
    % PATCHES VISUALISATION
    for lev=1:1:u_levels
        PLABELS = ULABELS(1:lev,:);
        [nLevels,nVertices] = size(PLABELS);
        PHI = HierarchicalInterface;
        PHI.nL = nLevels;
        
        % convert from clusterindex to listindex
        newLABELS = zeros(size(PLABELS));
        for l=1:size(PLABELS,1)
            for lm = 1:size(PLABELS,2)
                cl = PLABELS(l,lm);
                newLABELS(l,lm) = LC2Ind(PHI,l,cl);
            end
        end
        % For each point find the non-nan listindex
        VertexLabels = zeros(1,nVertices);
        for i=1:nVertices
            index = find(~isnan(PLABELS(:,i)));
            VertexLabels(i) = newLABELS(index(end),i);
        end
        [UV,~,VertexLabels] = unique(VertexLabels);
        
        RefScan = clone(preprocTemplate);
        RefScan.VertexValue = VertexLabels;
        RefScan.ColorMode = "Indexed";
        v = viewer(RefScan);
        v.Tag = 'All segments on Template';
        RefScan.ViewMode = "Solid";
        RefScan.Material = "Dull";
        colorbar(v.RenderAxes,'color',[1,1,1]);
        switch regionName
            case {'LH' 'SH'}
                view(-90,0)
            case 'RH'
                view(90,0)
        end
        v.SceneLightVisible = 1;
        v.SceneLightLinked = true;
        colormap(v.RenderAxes,'colorcube') % colorcube
        RefScan.PatchHandle.FaceColor = 'flat';
        print(v.Figure,'-dpng','-r300',[RESULTS_DIR  'segmentation_nL' num2str(nLevels)]);
    end
    clusterArray = LABELS;
else
    
    clustered = hierarchicalClustering(similarityMat,NUM_LEVELS,true,3,SEED);
    fig = paintClusters(clustered, preprocTemplate, NUM_LEVELS);
    savefig(fig, [SEGMENTATION_DIR  'segmentation.fig']);
    saveas(fig, [SEGMENTATION_DIR 'segmentation.png']);
    clusterArray = getClustersAsArray(clustered, NUM_LEVELS);
    save([SEGMENTATION_DIR 'segmentation.mat'],'clusterArray','-v7.3');
end

clear A similarityMat
%% Compute the PCA Components
% Investigate the number of PCA components to keep for each segment based on the explained Variance
% Get the partition of A that corresponds to the specific cluster with C
% points and then transform it into a N*3 x C array
partitions_num = 2^(NUM_LEVELS + 1) - 1;
explainedThresholds = 50:10:90;
clusterExpComponentsNum = zeros(partitions_num, length(explainedThresholds));
means = zeros(partitions_num, 3, nSubj);
templateCenters = zeros(partitions_num,3);

reshapedResT = reshape(resT', DIM, nVertices, nSubj);

parfor c=1:partitions_num
    try
        k = length(de2bi(c));
        clusterVec = clusterArray(k, :);
        i = c - (2^(k-1) - 1); % cluster index
        selectionMask = clusterVec==i;
        if ~any(selectionMask), continue; end
        ftResT = reshape(reshapedResT(:, selectionMask, :), DIM * sum(selectionMask), nSubj)';
        [~, pcaScores{c}, ~, ~, explained{c}] = pca(ftResT, 'Centered', true);
        explainedAccum = cumsum(explained{c});
        clusterExpComponentsNum(c,:) = sum(explainedAccum < explainedThresholds, 1);
        means(c, :, :) = mean(reshapedResT(:, selectionMask, :), 2);
        templateCenters(c, :) = mean(preprocTemplate.Vertices(selectionMask, :), 1);
    catch ME
        disp(c)
        rethrow(ME)
    end
end
fig = figure;
bar(explainedThresholds, clusterExpComponentsNum');
xlabel('Explained Variance Threshold');
ylabel('Number of Components Required');
title('Explained Variance~Number Of Components per Cluster');
savefig(fig, [SEGMENTATION_DIR 'explainedVariance.fig']);
saveas(fig, [SEGMENTATION_DIR 'explainedVariance.png']);
%%
selectedVarianceThreshold = 80;
maxNumComponents = max(clusterExpComponentsNum(: , explainedThresholds==selectedVarianceThreshold));
parfor k = 1:partitions_num
    scores = pcaScores{k};
    explainedAccum = cumsum(explained{k});
    clusterPCAPhenoFeatures{k} = scores(:, explainedAccum < selectedVarianceThreshold);
end
%%
save([SEGMENTATION_DIR, 'phenotype_varThres' num2str(selectedVarianceThreshold) '.mat'],'explained', 'clusterPCAPhenoFeatures','templateCenters','means','selectedVarianceThreshold','preprocPhenoIID', '-v7.3');


%%
% Show the behavior of the means of the clustered regions in 3D

localBehaviors = zeros(partitions_num,3, 3);
for i=1:partitions_num
    localBehaviors(i,:, :) = pca(squeeze(means(i, :, :))', 'NumComponents',3);
end
fig = figure;
ax = gca;
hold on
quiver3(ax, templateCenters(:,1), templateCenters(:,2), templateCenters(:,3), localBehaviors(:, 1, 1), localBehaviors(:, 2, 1), localBehaviors(:,3,1), 'LineWidth', 3);
quiver3(ax, templateCenters(:,1), templateCenters(:,2), templateCenters(:,3), localBehaviors(:, 1, 2), localBehaviors(:, 2, 2), localBehaviors(:,3,2), 'LineWidth', 3);
quiver3(ax, templateCenters(:,1), templateCenters(:,2), templateCenters(:,3), localBehaviors(:, 1, 3), localBehaviors(:, 2, 3), localBehaviors(:,3,3), 'LineWidth', 3);
legend('PC1', 'PC2', 'PC3');
hold off
brainShape = clone(preprocTemplate);
brainShape.RenderAxes = gca;
brainShape.Visible=true;
brainShape.ViewMode = 'SolidWireframe';
brainShape.Material = 'Dull';
axis equal
axis off
set(fig, 'color', 'none');
set(ax, 'color', 'none');
view(90, 0)
saveas(fig, [SEGMENTATION_DIR 'clusterPCA_interior.png']);
view(-90, 0)
saveas(fig, [SEGMENTATION_DIR 'clusterPCA_exterior.png']);
diary off;

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



