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
%
load([outpath 'SimilarityMatrix/SimilarityMatrixThickness']);

%% STEP 2: BUILDING RV MATRIX combine with distance map
[nVertices,DIM,nSubj] = size(Stage00ThicknessSym);
DistanceMatrix = nan*zeros(nVertices,nVertices);
for i=1:nVertices
    dist = intraDistancesv2(Region.AvgShape,'VertexIndex',i); % geodesic distance obj=Region.AvgShape;
    DistanceMatrix(i,:)= dist;
end
save([outpath '/SimilarityMatrix/DistanceMatrixAvgShape'],'DistanceMatrix','-v7.3'); 

DistanceMatrix_n = 1./DistanceMatrix;
DistanceMatrix_n(1:1+size(DistanceMatrix_n,1):end) = 1;  %diag = 1

%% STEP 3: RUNNING CLUSTERING
n_levels = 9;
type = 'weiss';
minPercValue = 1;
runs = 50;
% check how much percentage of DistanceMatrix
w = [0 0.1 0.2 0.25 0.3 0.55 0.8 1];
ConnectivityTab = zeros(length(w),4);
NMIscore = zeros(length(w),6);
for k=1:length(w)
    % weighted average 
    a = w(k);
    SimilarityMatrixThicknessDistMap = double(a*DistanceMatrix_n + (1-a)*SimilarityMatrixThickness);

    [LABELS,MASK] = HierarchicalFacialSegmentationv4(SimilarityMatrixThicknessDistMap,n_levels,type,runs,minPercValue);
    save([outpath 'Segmentation/WeissSegmentationThickness_' num2str(a*100) '%DistanceMap'],'LABELS','MASK','-v7.3');
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
    % PATCHES VISUALISATION
    figdir = [figpath 'Segmentation_Thickness_' num2str(a*100) '%DistanceMap'];
    mkdir(fullfile(figdir));

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
        print(v.Figure,'-dpng','-r300',[figdir '/Segmentation_Thickness_' num2str(a*100) '%DistanceMap_L' num2str(nLevels)]);
    end
    close all;
    %% METRICS
    % 1. Check connected component labeling 
    AdjacencyMatrix = full(Region.AvgShape.Adjacency);
    Connectivity = connectivityCheck(AdjacencyMatrix,LABELS,MASK);
    save([outpath '/Connectivity/ConnectivityThickness_' num2str(a*100) '%DistanceMap'],'Connectivity','-v7.3');
    ConnectivityTab(k,1) = a*100;
    ConnectivityTab(k,2) = length(find(VMASK));
    ConnectivityTab(k,3) = sum(Connectivity(:,2)==1);
    ConnectivityTab(k,4) = mean(Connectivity(:,3));
    % 2. NMIscore
    IN = load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/PHENOTYPING/SWH/Segmentation_nL9.mat');
    ULABELS_shape = IN.LABELS(index,:);
    for l=1:5 %i=1:u_levels
        disp(['Compute NMI for Level ' num2str(l)]);
        CLInd1 = ULABELS(l,:);
        CLInd2 = ULABELS_shape(l,:);
        [~,v] = nmi(CLInd1, CLInd2);
        NMIscore(k,1) = a*100;
        NMIscore(k,l+1) = v;
    end 
% % 3. Compare segement methods Silhouette score
% Sil= zeros(1,u_levels);
% for i= 2:u_levels
%     disp(['Compute Silhouette score for Level ' num2str(i)]);
%     CLInd2 = ULABELS2(i,:);
%     fn=@mydistfun;
%     s_list = silhouette(DistanceMatrix,CLInd2',fn);
%     %s = silhouette3d(DistanceMatrix, CLInd2);
%     Sil(i) = mean(s_list);
% end
end
ConnectivityTab_DistanceMap = array2table(round(ConnectivityTab,4));
ConnectivityTab_DistanceMap.Properties.VariableNames(1:4) = {'DistanceMapPct','TotalSegments','OnePieceSegments','MeanCoverage'};
writetable(ConnectivityTab_DistanceMap,[outpath '/ConnectivityTab_DistanceMap.xlsx']);  

NMITab_DistanceMap = array2table(round(NMIscore,4));
NMITab_DistanceMap.Properties.VariableNames(1:6) = {'DistanceMapPct','Level1','Level2','Level3','Level4','Level5'};
writetable(NMITab_DistanceMap,[outpath '/NMITab_DistanceMap.xlsx']);  
%%  END

%% visualize how a point correlates to the rest of the brain
point = 10000;
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = SimilarityMatrixT(point,:);
RefScan.ColorMode = "Indexed";
v1 = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
colorbar(v1.RenderAxes,'color',[1,1,1]);
v1.SceneLightVisible = 1;
v1.SceneLightLinked = true;
view(-90,0);
title(['Visualization of Thickness SimilarityMatrix at point ' num2str(point)]);
print(v1.Figure,'-dpng','-r300',['Corr_Thickness_' num2str(point)]);  
%
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = DistanceMatrix_n(point,:);
RefScan.ColorMode = "Indexed";
v2 = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
colorbar(v2.RenderAxes,'color',[1,1,1]);
v2.SceneLightVisible = 1;
v2.SceneLightLinked = true;
view(-90,0);
title(['Visualization of AvgShape 1/DistanceMatrix at point ' num2str(point)]);
print(v2.Figure,'-dpng','-r300',['Corr_1DistanceAvgShape_' num2str(point)]);  
%
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = SimilarityMatrixThicknessDistMap(point,:);
RefScan.ColorMode = "Indexed";
v5 = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
colorbar(v5.RenderAxes,'color',[1,1,1]);
v5.SceneLightVisible = 1;
v5.SceneLightLinked = true;
view(-90,0);
title(['Visualization of ThicknessSimilarityMatrix + ' newline num2str(a*100) '%AvgShape.DistanceMatrix at point ' num2str(point)]);
print(v5.Figure,'-dpng','-r300',['Corr_Thickness_' num2str(a*100) '%DistanceMap_' num2str(point)]); 
close all;
% visualize listindex = 8 ############################################################
% Found actually some shape segments are not one-piece
listindex = 8;
VertexLabels = zeros(1,nVertices);
[l,c] = Ind2LC(UHI,listindex);
CLInd = LABELS(l,:);
subInd = find(CLInd==c);% Find points that belong to listindex = 8
VertexLabels(subInd) = 1; 
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = VertexLabels;
RefScan.ColorMode = "Indexed";
v8 = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
%colorbar(v8.RenderAxes,'color',[1,1,1]);
v8.SceneLightVisible = 1;
v8.SceneLightLinked = true;
view(-90,0);
title(['Visualization of Segmentation based on Shape for ListIndex' num2str(listindex)]);
print(v8.Figure,'-dpng','-r300',['Segmentation_ListIndex' num2str(listindex)]); 

