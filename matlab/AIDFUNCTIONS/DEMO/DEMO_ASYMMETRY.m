%% Investigating LEFT - RIGHT asymmetry
close all;clear;
%%
% DATA_DIR = '../SAMPLE_DATA/';
% THREADS = 8;
% nSamples = 100;
% reduce = 0.01;
% nRep = 3;
% performExperiments=1;
% nIter = 1000;
DATA_DIR = '/';
THREADS = 40;
nSamples = 1000;
reduce = 0.01;
nRep = 3;
nIter = 10000;

performExperiments = 0;
restoredefaultpath;

addpath(genpath('AIDFUNCTIONS'));
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
Render = cell(1,2);
for r=1:nR
    %r=2
    disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
    regphenopath = [phenopath Regions{r} '/'];
    DATA{r} = load([regphenopath 'STAGE00DATA']);
    Render{r} = load([regphenopath 'RENDERMATERIAL.mat']);
    
end

%% Subsampling space to match R version which has memory issues to pick the whole space

Ns = floor(1/reduce);
pdim = size(DATA{1}.Region.AlignedShapes,1);
LH = DATA{1}.Region.AlignedShapes(1:Ns:end, :, 1:nSamples);
RH = DATA{2}.Region.AlignedShapes(1:Ns:end, :, 1:nSamples);
RH(:,1,:,:) = -1*RH(:,1,:,:);
Template = clone(DATA{1}.Region.AvgShape);
% The folowing will not work as the number of vertices differs per reduced
% template
% if reduce ~= 1
%     assert(reduce<1)
%     facesnum =  floor(size(LH,1)*reduce) * 2; % the faces number is ~ 2 * vertices
%     faces = Template.Faces;
%     reducedTemplateStruct = reducepatch(...
%         struct("vertices", Template.Vertices, "faces", Template.Faces) ,  facesnum);    
%     vernum = size(reducedTemplateStruct.vertices,1);
%     reducedLH = zeros(vernum,size(LH,2),size(LH,3));
%     for i=1:size(reducedLH,3)
%         reducedLH(:,:,i) = reducepatch(...
%         struct("vertices", LH(:,:,i), "faces", faces) ,  facesnum).vertices;
%     end
% 
% else
%     reducedTemplateStruct.vertices = Template.Vertices;
%     reducedTemplateStruct.faces = Template.Faces;
% end


reducedTemplate = shape3D;
reducedTemplate.struc2obj(Template.obj2struc());

reducedTemplate.Vertices = reducedTemplate.Vertices(1:Ns:end,:);
reducedTemplateIndices = 1:Ns:pdim;
reducedTemplateAdjacency = Template.Adjacency;

reducedTemplate.Vertices = reducedTemplate.Vertices - mean(reducedTemplate.Vertices, 1);
reducedTemplate.Vertices = (reducedTemplate.Vertices ./ sqrt(sum(reducedTemplate.Vertices.^2,'all')));
strReducedTemplate = reducedTemplate.obj2struc();
strTemplate = Template.obj2struc();
%% Prepare Data to be provided to GPA

% figure;
% hold on;
% scatter3(RH(:,1,1),RH(:,2,1),RH(:,3,1),'b');,
% scatter3(LH(:,1,1),LH(:,2,1),LH(:,3,1),'r');,
% hold off;


%% prepare output to R

% 
% 
% nRep = 6;
% % 
% input_to_r_path =  [DATA_DIR 'm2r.mat'];
% save(input_to_r_path, 'LH', 'RH','strReducedTemplate', 'strTemplate', 'nSamples','Ns','nRep');

% TotalRepShapes: 

%% GPA
% input_from_r_path = [DATA_DIR 'r2m.mat'];
% gpa_from_r = load(input_from_r_path);
%%
 
TotalShapes = cat(3,LH,RH);
%%
% [AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(TotalShapes,Template,3,true,false,true,false);
[AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(TotalShapes,reducedTemplate,3,true,'best',true,false);





%%
% clear DATA;
clear  LH RH TotalShapes;
%% TWO WAY PROCRUSTES ANOVA ON REDUCED DATA
LHAligned = AlignedShapes(:,:,1:size(AlignedShapes, 3)/2 );
RHAligned = AlignedShapes(:,:,size(AlignedShapes, 3)/2+1:end);
%%
%%
% LHAligned_R=gpa_from_r.LHAligned;
% RHAligned_R = gpa_from_r.RHAligned;
% reducedTemplateAligned_R = gpa_from_r.TemplateAligned;
%%
% 
% display3DLandmarks(reducedTemplate, reducedTemplateAligned_R);
% %%
% display3DLandmarks(LHAligned(:,:,1), LHAligned_R(:,:,1));
% 
% %%
% display3DLandmarks(reducedTemplateAligned_R, LHAligned_R(:,:,1));
% %%
% 
% display3DLandmarks(LHAligned(:,:,1), reducedTemplate.Vertices);
%%


% figure;
% hold on;
% scatter3(RHAligned(:,1,1),RHAligned(:,2,1),RHAligned(:,3,1),'b');,
% scatter3(LHAligned (:,1,1),LHAligned (:,2,1),LHAligned (:,3,1),'r');,
% hold off;


% LHAlignedInt16 = int16(LHAligned.*10000);clear LHAligned;
% RHAlignedInt16 = int16(RHAligned.*10000);clear RHAligned;


Shapes = permute(AlignedShapes,[2 1 3]);
Shapes = reshape(Shapes,size(Shapes,1)*size(Shapes,2),size(Shapes,3))';
clear AlignedShapes;

%%
RepShapes = zeros(size(Shapes,1),size(Shapes,2),nRep,'single');
mag = var(Shapes,0,2);
for i=1:1:nRep
    RepShapes(:,:,i) = single(Shapes) + single(randn(size(Shapes,1),size(Shapes,2)).*mag*0.2);
end
mult = double(intmax('int16')) / (max(abs(RepShapes),[],'all'));
RepShapesInt16 = int16(RepShapes.*mult);

figure; histogram(reshape(RepShapesInt16 - int16(mult * Shapes), 1,[])); title({'Landmark Coordinate dislocation','for generated replications'});
%%
clear RepShapes;
%%
X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:); 
%%
out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,nIter, mult);

%%

if performExperiments
    toCheckSizes= [10, 20, 50, 100];
    retS = checkSize(X1, X2, toCheckSizes);
    plotExp(retS, toCheckSizes, 'Population Size');

    %%
    toCheckLSizes= [50, 100, 150, 200, 250];
    retL = checkLSize(X1, X2, toCheckLSizes);
    plotExp(retL, toCheckLSizes, 'Number of Landmarks');

    %%

    toCheckRepsNum= [1,2,3];
    retR = checkRepsNum(X1, X2, toCheckRepsNum);
    plotExp(retR, toCheckRepsNum, 'Number of Replications');

    %%
    toCheckPermsNum= [50,100,180,300];
    retP = checkPermsNum(X1, X2, toCheckPermsNum);
    plotExp(retP, toCheckPermsNum, 'Number of Permutations');
end
%% Upsampling
toupsample = [out.LM.I;out.LM.IF;out.LM.permIF;out.LM.D;out.LM.DF; out.LM.permDF; out.LM.F;out.LM.FF;out.LM.permFF]';

outupsampled= upsampleShape3D(toupsample, reducedTemplateAdjacency, reducedTemplateIndices);
outu.LM.I = outupsampled(:,1);
outu.LM.IF = outupsampled(:,2);
outu.LM.permIF = outupsampled(:,3);
outu.LM.D = outupsampled(:,4);
outu.LM.DF = outupsampled(:,5);
outu.LM.permDF = outupsampled(:,6);
outu.LM.F = outupsampled(:,7);
outu.LM.FF= outupsampled(:,8);
outu.LM.permFF = outupsampled(:,9);
%%
showstruct = outu;
 %% BELOW IS AN IDEA OF RENDERING, BUT WILL NOT WORK BECAUSE WE DO NOT HAVE ALL THE MESH POINTS

i=1;
VertexValues{i,1} = showstruct.LM.I;titlenames{i,1} = 'I';
VertexValues{i,2} = showstruct.LM.IF;titlenames{i,2} = 'IF';
VertexValues{i,3} = showstruct.LM.permIF;titlenames{i,3} = 'p-Value';
i=i+1;
VertexValues{i,1} = showstruct.LM.D;titlenames{i,1} = 'D';
VertexValues{i,2} = showstruct.LM.DF;titlenames{i,2} = 'DF';
VertexValues{i,3} = showstruct.LM.permDF;titlenames{i,3} = 'p-Value';
i=i+1;
VertexValues{i,1} = showstruct.LM.F;titlenames{i,1} = 'F';
VertexValues{i,2} = showstruct.LM.FF;titlenames{i,2} = 'FF';
VertexValues{i,3} = showstruct.LM.permFF;titlenames{i,3} = 'p-Value';
labels = ["Individual", "Directional", "Fluctuating"];
thresholds = [0.05, 0.01, 0.005, 0.0001];
totalStats = out.Total;
brainSurface = Render{1};
nValues = 3;
for i=1:nValues
    val = VertexValues{i,3};
    res = zeros(size(val));
    c = 1;
    for t=thresholds 
        res(val<=t) = c;
        c = c+1;
    end
    VertexValues{i,4} = res;
    titlenames{i,4} = "Significance";
end
data = VertexValues;
save('asymmetry.mat', "data", "titlenames","labels","thresholds","totalStats","brainSurface",'-v7.3');
%%
load('asymmetry.mat');
f = figure;f.Position = [95  98  2192  1106];f.Color = [1 1 1];%

arrange = [3 8];
counter = 0;
clim = [];
rend = brainSurface;
cmap = parula(256);


dmap = parula(length(thresholds));


thresholdsTicks = thresholds;
for i=1:nValues
    for j=1:4
        if j==4, map=dmap; else, map=cmap; end
        %i=1;
       counter = counter+1;
       fout.ax1{i,j} = subplot(arrange(1),arrange(2),counter,'Parent',f);
       colormap(fout.ax1{i,j},map);
       
       if ~isempty(clim), set(fout.ax1{i,j},'clim',clim);end
       renderBrainSurface(rend,VertexValues{i,j},fout.ax1{i,j});
       colorbar(fout.ax1{i,j},'SouthOutside');
       if j<4        
           set(fout.ax1{i,j},'clim',[0 max(VertexValues{i,j})]);
       end
       if j==4,set(fout.ax1{i,j},'clim',[0 length(thresholds)]); cb=findall(gcf,'type','ColorBar');cb(1).Ticks=1:length(thresholdsTicks);
           cb(1).TickLabels=thresholdsTicks; cb(1).Label.String = 'p-value'; end
       view(fout.ax1{i,j},rend.viewval(1),0);
       light = camlight(fout.ax1{i,j},'headlight');
       set(light,'Position',get(fout.ax1{i,j},'CameraPosition'));
       drawnow;
       title(fout.ax1{i,j},titlenames{i,j})
       counter = counter+1;
       fout.ax2{i,j} = subplot(arrange(1),arrange(2),counter,'Parent',f);
       renderBrainSurface(rend,VertexValues{i,j},fout.ax2{i,j});
       view(fout.ax2{i,j},-1*rend.viewval(1),0);
       colorbar(fout.ax2{i,j},'SouthOutside');
       if j<4
           set(fout.ax2{i,j},'clim',[0 max(VertexValues{i,j})]);
       end
       colormap(fout.ax2{i,j},map);
       if j==4,set(fout.ax2{i,j},'clim',[0 length(thresholds)]); cb=findall(gcf,'type','ColorBar');cb(1).Ticks=1:length(thresholdsTicks);
           cb(1).TickLabels=thresholdsTicks; cb(1).Label.String = 'p-value'; end
       light = camlight(fout.ax2{i,j},'headlight');
       set(light,'Position',get(fout.ax2{i,j},'CameraPosition'));
       drawnow;
       if ~isempty(clim), set(fout.ax2{i,j},'clim',clim);end
    end
end
saveas(f, ['../results/demo_asymmetry/results_' num2str(nSamples) '_' num2str(Ns) '.fig']);
%%
% load('asymmetry.mat');
% for i=1:4
%         map=dmap;
%        counter = counter+1;
%        viewer = viewer3D;
%        viewer.set("Tag", char(labels(i)));
%        shape = rend.RefScan;
%        shape.TextureMap =  dmap;
%        shape.ColorMode = "indexed";
%        shape.ViewMode = "solid";
%        shape.VertexValue = VertexValues{i,4};
%        shape.viewer(viewer);
% end
