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

%% Smooth thickness 
% /IMAGEN/AIDFUNCTIONS/ObjectClasses/3D Objects/@meshObj/...
% edit curvature
% edit signedCurvature   %Patch Curvature by Dirk-Jan Kroon
% edit smoothFunction
% gradient = obj.VertexNormals;
w = [2 5 10 20 30 50 100];
ConnectivityTab = zeros(length(w),4);
NMIscore = zeros(length(w),6);
for k=1:length(w) 
    iter = w(k);
    SmoothThickness2DMatrix = nan*zeros(nVertices,nSubj);
    for s=1:nSubj
        thickness = double(Thickness2DMatrix(:,s));
        obj = clone(Region.AvgShape);
        obj.Vertices = Region.AlignedShapes(:,:,s);
        smoothThickness = smoothFunctionv2(obj,thickness,iter,'functiondistance'); 
        SmoothThickness2DMatrix(:,s) = single(smoothThickness);
    end   
    avgT = mean(SmoothThickness2DMatrix,2);
    resT = getResiduals([COV.DATA Region.CentroidSizes(:)],SmoothThickness2DMatrix');
    resT = repmat(avgT',size(resT,1),1)+resT;
    save([outpath '/CorrectedData/COVThicknessSmooth_iter' num2str(iter)],'resT','-v7.3');
    %% STEP 2: BUILDING RV MATRIX
    type = 'cov';
    SimilarityMatrixSmoothThickness = buildRVmatrixDim(resT,type,1);
    save([outpath '/SimilarityMatrix/SimilarityMatrixSmoothThickness_iter' num2str(iter)],'SimilarityMatrixSmoothThickness','-v7.3'); 
    %% STEP 3: RUNNING CLUSTERING
    n_levels = 9;
    type = 'weiss';
    minPercValue = 1;
    runs = 50;
    [LABELS,MASK] = HierarchicalFacialSegmentationv4(SimilarityMatrixSmoothThickness,n_levels,type,runs,minPercValue);
    save([outpath 'Segmentation/WeissSegmentationThickness_Smooth_iter' num2str(iter)],'LABELS','MASK','-v7.3');
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
    load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/PHENOTYPING/SWH/STAGE00DATA.mat');
    figdir = [figpath 'Segmentation_Thickness_Smooth_iter' num2str(iter)];
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
        print(v.Figure,'-dpng','-r300',[figdir '/Segmentation_Thickness_Smooth_iter' num2str(iter) '_L' num2str(nLevels)]);
    end
    close all;
    %% METRICS
    % 1. Check connected component labeling 
    AdjacencyMatrix = full(Region.AvgShape.Adjacency);
    Connectivity = connectivityCheck(AdjacencyMatrix,LABELS,MASK);
    save([outpath '/Connectivity/ConnectivityThickness_Smooth_iter' num2str(iter)],'Connectivity','-v7.3');
    ConnectivityTab(k,1) = iter;
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
        NMIscore(k,1) = iter;
        NMIscore(k,l+1) = v;
    end 
end 
ConnectivityTab_Smooth = array2table(round(ConnectivityTab,4));
ConnectivityTab_Smooth.Properties.VariableNames(1:4) = {'SmoothIter','TotalSegments','OnePieceSegments','MeanCoverage'};
writetable(ConnectivityTab_Smooth,[outpath '/ConnectivityTab_Smooth.xlsx']);  

NMITab_Smooth = array2table(round(NMIscore,4));
NMITab_Smooth.Properties.VariableNames(1:6) = {'SmoothIter','Level1','Level2','Level3','Level4','Level5'};
writetable(NMITab_Smooth,[outpath '/NMITab_Smooth.xlsx']);  
%%  END

%% STUDY 
obj = Region.AvgShape;
FV.vertices = obj.Vertices;
FV.faces = obj.Faces;
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(FV);
% visualize gaussian/mean curvature
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = Cgaussian; %Cmean
RefScan.ColorMode = "Indexed";
v3 = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
colorbar(v3.RenderAxes,'color',[1,1,1]);
v3.SceneLightVisible = 1;
v3.SceneLightLinked = true;
view(-90,0);
title('Visualization of gaussian curvature');
print(v3.Figure,'-dpng','-r300','GaussianCurvature_AvgShape');  
%% Smooth curve
sCurvature = signedCurvaturev2(obj);
% visualize Smoothed curvature
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = sCurvature;
RefScan.ColorMode = "Indexed";
v = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
colorbar(v.RenderAxes,'color',[1,1,1]);
v.SceneLightVisible = 1;
v.SceneLightLinked = true;
view(-90,0);
title('Visualization of smoothed curvature');
print(v.Figure,'-dpng','-r300','SmoothedCurvature_AvgShape'); 
%% Smoothing of thickness 
avgT = double(mean(Thickness2DMatrix,2));
% visualize average thickness
RefScan = clone(Region.AvgShape);
RefScan.VertexValue = avgT;
RefScan.ColorMode = "Indexed";
v1 = viewer(RefScan);
RefScan.ViewMode = "Solid";
RefScan.Material = "Dull";
colorbar(v1.RenderAxes,'color',[1,1,1],'Limits',[0.02 0.14]);
v1.SceneLightVisible = 1;
v1.SceneLightLinked = true;
view(-90,0);
title('Visualization of average thickness');
print(v1.Figure,'-dpng','-r300',[figpath 'SmoothThickness/AvgThicknessPlotOnAvgShape']); 
% check iteration
w = [2 5 10 20 30 50 100];
for k=1:length(w) 
    iter = w(k);
    obj = clone(Region.AvgShape);
    smoothThickness = smoothFunctionv2(obj,avgT,iter,'functiondistance'); 
    % visualize Smoothed thickness
    RefScan = clone(Region.AvgShape);
    RefScan.VertexValue = smoothThickness;
    RefScan.ColorMode = "Indexed";
    v2 = viewer(RefScan);
    RefScan.ViewMode = "Solid";
    RefScan.Material = "Dull";
    colorbar(v2.RenderAxes,'color',[1,1,1],'Limits',[0.02 0.14]);  
    v2.SceneLightVisible = 1;
    v2.SceneLightLinked = true;
    view(-90,0);
    title(['Visualization of smoothed average thickness at iter = ' num2str(iter)]);
    print(v2.Figure,'-dpng','-r300',[figpath 'SmoothThickness/SmoothedThicknessPlotOnAvgShape_iter' num2str(iter)]); 
end
close all;







