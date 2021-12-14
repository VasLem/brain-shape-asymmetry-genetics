clc;clear all;
addpath(genpath('/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/SHARED/AIDFUNCTIONS/'));
addpath(genpath('/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/FUNCTIONS/'));
 
outpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/01_Pheno/Data/';
figpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/01_Pheno/Fig/';
%load average thickness L/R
load([outpath '/Stage00Data/Stage00ThicknessSym.mat']);
%load COV and SHAPE
load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/COVARIATES/COVDATAINLIERS.mat');
load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/PHENOTYPING/SWH/STAGE00DATA.mat');

%% STEP 1: CORRECTING THE DATA, Correction for covariate Data % + BrainSize
[nVertices,DIM,nSubj] = size(Stage00ThicknessSym);
Thickness2DMatrix = permute(Stage00ThicknessSym,[2 1 3]);
Thickness2DMatrix = reshape(Thickness2DMatrix,nVertices*DIM,nSubj);
Thickness2DMatrix = Thickness2DMatrix./Region.CentroidSizes(:)'; % scale thickness

avgT = mean(Thickness2DMatrix,2);
resT = getResiduals([COV.DATA Region.CentroidSizes(:)],Thickness2DMatrix');
resT = repmat(avgT',size(resT,1),1)+resT;
save([outpath '/CorrectedData/COVThickness'],'resT','-v7.3');

%% STEP 2: BUILDING RV MATRIX
type = 'cov';
SimilarityMatrixThickness = buildRVmatrixDim(resT,type,1);
save([outpath 'SimilarityMatrix/SimilarityMatrixThickness'],'SimilarityMatrixThickness','-v7.3'); 

%% STEP 3: RUNNING CLUSTERING
n_levels = 9;
type = 'weiss';%'symmetric laplacian'; %'ratiocute';%'ncute';%'weiss'; is not so 'symmetric laplacian'
    %nice for RH and SH;
    %type = 'weiss';%'ratiocute';%'ncute';%'weiss';
minPercValue = 1;
runs = 50;
[LABELS,MASK] = HierarchicalFacialSegmentationv4(SimilarityMatrixThickness,n_levels,type,runs,minPercValue);
save([outpath 'Segmentation/WeissSegmentationThickness'],'LABELS','MASK','-v7.3');

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
figdir = [figpath 'Segmentation_Thickness'];
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
        print(v.Figure,'-dpng','-r300',[figdir '/Segmentation_Thickness_nL' num2str(nLevels)]);
    end
close all;
%% END

%% Average thickness for each segement
% For each cluster find the consisting Vertices
clusterInd = find(Render.UMASK);
cluster_vertex = zeros(length(clusterInd), nVertices);
clusterAvgThickness = zeros(length(clusterInd),1);
for cI=1:length(clusterInd)
    [~, col]= find(newLABELS == clusterInd(cI));
    cluster_vertex(cI,col)= 1;
    clusterAvgThickness(cI)= sum(AvgThickness_S(col))/length(col);
end

