%% Investigating LEFT - RIGHT asymmetry
close all;clear;
restoredefaultpath;
%%
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('FUNCTIONS'));

%%
loadWhileInLab = 1;

seed = 42;
rng(seed); % For reproducible results
if endsWith(cd, "AIDFUNCTIONS/DEMO")
    cd("../..")
end
DATA_DIR = '../SAMPLE_DATA/';
RESULTS_DIR = '../results/hierarchicalClusteringDemo/';
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
%%
applyOnAtlas = false;
% THREADS = 8;atlas
nPicks = 5;
nSamplesPerPick = 50;
% nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
% reduce = 0.05;
nRep = 3;
nIter = 1000;
THREADS = 32;
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
covGenoPath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/COVARIATES/COVDATA2USE.mat'];

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
Region =  DATA{1}.Region;
%% 
brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
refTemplate = brainSurface.RefScan;



%%
[preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, 1, subsample, 3);

nLandmarks = length(preprocLandmarksIndices);
%%
selection = 'Asymmetry';
switch selection
    case 'Asymmetry'
        A = (preprocLH - preprocRH);
    case 'Symmetry'
        A =  (preprocLH + preprocRH)/2;
end

%% STEP 1:
[nVertices,DIM,nSubj] = size(A);
sym2DMatrix = permute(A,[2 1 3]);
sym2DMatrix = reshape(sym2DMatrix,nVertices*DIM,nSubj);
disp("Applying covariates effect subtraction (control) from phenotype..")
% Align covariates with phenotype
covariates = load(covGenoPath).COV;
%%
covAssignmentMatrix = (str2double(phenoIID) == str2double(covariates.IID)');
[covPhenoIndex, covIndex] = find(covAssignmentMatrix);

covData = covariates.DATA(covIndex, :);
sym2DMatrix = sym2DMatrix(:, covPhenoIndex);
phenoIID = phenoIID(covPhenoIndex);

avgT = mean(sym2DMatrix,2);
% Removing components that can be controlled from covariates, keeping only
% the residuals from the corresponding PLS model
%%
disp("Fitting PLS model to covariates and removing their effect from the phenotype..");
resT = getResiduals(covData,sym2DMatrix');
resT = repmat(avgT',size(resT,1),1)+resT;
disp("Saving computed residuals..");
save([RESULTS_DIR 'COV' selection],'resT','-v7.3');
%% STEP 2: BUILDING RV MATRIX
disp("Computing RV matrix..")
similarityMat = buildRVmatrixDim(resT, 'cov', 3);
disp("Saving RV Matrix..")
save([RESULTS_DIR 'COVSimilarityMatrix' selection],'similarityMat','-v7.3');

% load([RESULTS_DIR 'Segmentation/SimilarityMatrix' selection]);
%% STEP 3: RUNNING CLUSTERING
n_levels = 9;
type = 'weiss';%'symmetric laplacian'; %'ratiocute';%'ncute';%'weiss'; is not so 'symmetric laplacian'
    %nice for RH and SH;
    %type = 'weiss';%'ratiocute';%'ncute';%'weiss';
minPercValue = 1;
runs = 50;
[LABELS,MASK] = HierarchicalFacialSegmentationv4(double(similarityMat),n_levels,type,runs,minPercValue);
save([RESULTS_DIR 'COVWeissSegmentation' selection],'LABELS','MASK','-v7.3');

%% VISUALIZING THE SEGMENTATION
% find the deepest still meaningfull level
levels = nansum(LABELS,2);
index = find(levels);
u_levels = index(end);% levels to analyze
ULABELS = LABELS(index,:);
UHI = HierarchicalInterface;
UHI.nL = u_levels;
UMASK = MASK(1:UHI.nLC);
disp(['Number of total clusters: ' num2str(length(find(UMASK)))]);
v_levels = 9;% ter visualisatie
VLABELS = ULABELS(1:v_levels,:);
VHI = HierarchicalInterface;
VHI.nL = v_levels;
VMASK = UMASK(1:VHI.nLC);
disp(['Number of visualized clusters: ' num2str(length(find(VMASK)))]);  

% make figdir
figdir = [RESULTS_DIR 'Segmentation' selection];
mkdir(fullfile(figdir));    
 
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

        RefScan = clone(Region.AvgShape);
        RefScan.VertexValue = VertexLabels;
        RefScan.ColorMode = "Indexed";
        v = viewer(RefScan);
        v.Tag = 'All segments on Template';
        RefScan.ViewMode = "Solid";
        RefScan.Material = "Dull";
        colorbar(v.RenderAxes,'color',[1,1,1]);
        switch Region.Name
            case {'LH' 'SH'} 
                view(-90,0)
            case 'RH'
                view(90,0)
        end
        v.SceneLightVisible = 1;
        v.SceneLightLinked = true;
        colormap(v.RenderAxes,'colorcube') % colorcube
        RefScan.PatchHandle.FaceColor = 'flat';
        print(v.Figure,'-dpng','-r300',[figdir '/Segmentation_' selection '_nL' num2str(nLevels)]);
    end

%%SAMPLE_DATA
% clear DATA;
%%

numLevels = 4;
clustered = hierarchicalClustering(similarityMat,numLevels,true,3,seed);
%%
fig = paintClusters(clustered, preprocTemplate, numLevels);
savefig(fig, [RESULTS_DIR 'summaryMySegmentation' selection '.fig']);
saveas(fig, [RESULTS_DIR 'summaryMySegmentation' selection '.png']);

%%
clusterVec = getClustersAsVec(clustered);
%%
clusters = unique(clusterVec);


%% Investigate the number of PCA components to keep for each segment based on the explained Variance
% Get the partition of A that corresponds to the specific cluster with C
% points and then transform it into a N*3 x C array

explainedThresholds = 80:4:96;
clusterExpComponentsNum = zeros(length(clusters), length(explainedThresholds));
parfor i=1:length(clusters)
    selectionMask = clusterVec==clusters(i);
    ftA = reshape(permute(A(selectionMask, :, :), [2, 1, 3]), size(A,2) * sum(selectionMask), size(A, 3))';
    [~, scores, ~, ~, explained] = pca(ftA, 'Centered', true);
    explainedAccum = cumsum(explained);
    clusterExpComponentsNum(i,:) = sum(explainedAccum < explainedThresholds, 1);
end
fig = figure;
bar(explainedThresholds, clusterExpComponentsNum');
xlabel('Explained Variance Threshold');
ylabel('Number of Components Required');
title('Explained Variance~Number Of Components per Cluster');
savefig(fig, [RESULTS_DIR 'explainedVariance_' selection '.fig']);
saveas(fig, [RESULTS_DIR 'explainedVariance_' selection '.png']);
%%
selectedVarianceThreshold = 88;
maxNumComponents = max(clusterExpComponentsNum(: , explainedThresholds==selectedVarianceThreshold));
means = zeros(length(clusters), 3, nSamples);
templateCenters = zeros(length(clusters),3);
parfor i=1:length(clusters)
    selectionMask = clusterVec==clusters(i);
    means(i, :, :) = mean(A(clusterVec == clusters(i), :, :), 1);
    templateCenters(i, :) = mean(template.Vertices(clusterVec==clusters(i), :), 1);
    ftA = reshape(permute(A(selectionMask, :, :), [2, 1, 3]), size(A,2) * sum(selectionMask), size(A, 3))';
    [~, scores, ~, ~, explained] = pca(ftA, 'Centered', true, 'NumComponents', maxNumComponents);
    explainedAccum = cumsum(explained);
    clusterPCAPhenoFeatures{i} = scores(:, explainedAccum < selectedVarianceThreshold);
end
%%
save([RESULTS_DIR, selection 'Phenotype.mat'], 'clusterVec','clusterPCAPhenoFeatures','templateCenters','means','selectedVarianceThreshold','phenoIID', '-v7.3');


%%
% Show the behavior of the means of the clustered regions in 3D

localBehaviors = zeros(length(clusters),3, 3);
for i=1:length(clusters)
    localBehaviors(i,:, :) = pca(squeeze(means(clusters(i), :, :))', 'NumComponents',3);
end
fig = figure;
ax = gca;
hold on
quiver3(ax, templateCenters(:,1), templateCenters(:,2), templateCenters(:,3), localBehaviors(:, 1, 1), localBehaviors(:, 2, 1), localBehaviors(:,3,1), 'LineWidth', 3);
quiver3(ax, templateCenters(:,1), templateCenters(:,2), templateCenters(:,3), localBehaviors(:, 1, 2), localBehaviors(:, 2, 2), localBehaviors(:,3,2), 'LineWidth', 3);
quiver3(ax, templateCenters(:,1), templateCenters(:,2), templateCenters(:,3), localBehaviors(:, 1, 3), localBehaviors(:, 2, 3), localBehaviors(:,3,3), 'LineWidth', 3);
legend('PC1', 'PC2', 'PC3');
hold off
brainShape = clone(template);
brainShape.RenderAxes = gca;
brainShape.Visible=true;
brainShape.ViewMode = 'SolidWireframe';
brainShape.Material = 'Dull';
axis equal
axis off
set(fig, 'color', 'none');
set(ax, 'color', 'none');
view(90, 0)
saveas(fig, [RESULTS_DIR 'clusterPCA_' selection '_interior.png']);
view(-90, 0)
saveas(fig, [RESULTS_DIR 'clusterPCA_' selection '_exterior.png']);


%%

function clusterVec = getClustersAsVec(clustered)
    clusterVec = ones(1,length(clustered.indices));
    if ~isempty(clustered.parts)
        part1 = getClustersAsVec(clustered.parts{1});
        part2 = getClustersAsVec(clustered.parts{2}) + max(part1);
        clusterVec(clustered.parts{1}.indices) = part1;
        clusterVec(clustered.parts{2}.indices) = part2;
    end
end



