clc;clear all;
addpath(genpath('/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/SHARED/AIDFUNCTIONS/'));
addpath(genpath('/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/FUNCTIONS/'));
 
outpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/01_Pheno/Data/';
figpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/01_Pheno/Fig/';

%
load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/PHENOTYPING/SWH/STAGE00DATA.mat');
load([outpath 'SimilarityMatrix/SimilarityMatrixThickness']);
load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/PHENOTYPING/SWH/SimilarityMatrix.mat');
%% STEP 3: RUNNING CLUSTERING
n_levels = 9;
type = 'weiss';
minPercValue = 1;
runs = 50;
% check how much percentage of DistanceMatrix
w = [0 0.1 0.2 0.3 0.4 0.5 0.8 1];
ConnectivityTab = zeros(length(w),4);
NMIscore = zeros(length(w),6);
for k=1:length(w)
    % weighted average 
    a = w(k);
    SimilarityMatrixThicknessShapeMap = double(a*SimilarityMatrix + (1-a)*SimilarityMatrixThickness);

    [LABELS,MASK] = HierarchicalFacialSegmentationv4(SimilarityMatrixThicknessShapeMap,n_levels,type,runs,minPercValue);
    save([outpath 'Segmentation/WeissSegmentationThickness_' num2str(a*100) '%ShapeMap'],'LABELS','MASK','-v7.3');
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
    figdir = [figpath 'Segmentation_Thickness_' num2str(a*100) '%ShapeMap'];
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
        print(v.Figure,'-dpng','-r300',[figdir '/Segmentation_Thickness_' num2str(a*100) '%ShapeMap_L' num2str(nLevels)]);
    end
    close all;
    %% METRICS
    % 1. Check connected component labeling 
    AdjacencyMatrix = full(Region.AvgShape.Adjacency);
    Connectivity = connectivityCheck(AdjacencyMatrix,LABELS,MASK);
    save([outpath '/Connectivity/ConnectivityThickness_' num2str(a*100) '%ShapeMap'],'Connectivity','-v7.3');
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
end
ConnectivityTab_ShapeMap = array2table(round(ConnectivityTab,4));
ConnectivityTab_ShapeMap.Properties.VariableNames(1:4) = {'ShapeMapPct','TotalSegments','OnePieceSegments','MeanCoverage'};
writetable(ConnectivityTab_ShapeMap,[outpath '/ConnectivityTab_ShapeMap.xlsx']);  

NMITab_ShapeMap = array2table(round(NMIscore,4));
NMITab_ShapeMap.Properties.VariableNames(1:6) = {'ShapeMapPct','Level1','Level2','Level3','Level4','Level5'};
writetable(NMITab_ShapeMap,[outpath '/NMITab_ShapeMap.xlsx']);  
%%  END
% Check connected component labeling 
AdjacencyMatrix = full(Region.AvgShape.Adjacency);
k = find(AdjacencyMatrix(1,:));
point = 100;
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = AdjacencyMatrix(point,:);
RefScan.ColorMode = "Indexed";
v4 = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
colorbar(v4.RenderAxes,'color',[1,1,1]);
v4.SceneLightVisible = 1;
v4.SceneLightLinked = true;
