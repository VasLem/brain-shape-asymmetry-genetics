%% Investigating LEFT - RIGHT asymmetry
close all;clear;
restoredefaultpath;
% addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));

addpath(genpath('AIDFUNCTIONS'));

%% SETTING UP COMPUTATION POWER
try
%     parpool('LocalSingle',20);
    parpool('local',6);
catch
end
%% GETTING SOME INFO ON THE BRAIN TEMPLATE
% DATA_DIR = '/';
DATA_DIR = '../SAMPLE_DATA/';

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
% nSamples = 1000;

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

nSamples = 100;
reduce = 0.01;
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
input_from_r_path = [DATA_DIR 'r2m.mat'];
gpa_from_r = load(input_from_r_path);
%%
 
TotalShapes = cat(3,LH,RH);
%%
% [AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(TotalShapes,Template,3,true,false,true,false);
[AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(TotalShapes,reducedTemplate,3,true,'best',true,false);





%%
clear DATA LH RH TotalShapes;
%% TWO WAY PROCRUSTES ANOVA ON REDUCED DATA
LHAligned = AlignedShapes(:,:,1:size(AlignedShapes, 3)/2 );
RHAligned = AlignedShapes(:,:,size(AlignedShapes, 3)/2+1:end);
%%
%%
LHAligned_R=gpa_from_r.LHAligned;
RHAligned_R = gpa_from_r.RHAligned;
reducedTemplateAligned_R = gpa_from_r.TemplateAligned;
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


LHAlignedInt16 = int16(LHAligned.*10000);clear LHAligned;
RHAlignedInt16 = int16(RHAligned.*10000);clear RHAligned;


Shapes = permute(AlignedShapes,[2 1 3]);
Shapes = reshape(Shapes,size(Shapes,1)*size(Shapes,2),size(Shapes,3))';
clear AlignedShapes;

%%
nRep = 6;
% nRep = 3;
RepShapes = zeros(size(Shapes,1),size(Shapes,2),nRep,'single');
var(Shapes);
for i=1:1:nRep
    RepShapes(:,:,i) = single(Shapes) + single(randn(size(Shapes,1),size(Shapes,2)).*var(Shapes,0,2)*0.5);
end
RepShapesInt16 = int16(RepShapes.*10000);clear RepShapes;
%%
X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:); 
%%
out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,1000);
%%
shape = shape3D;
shape.Vertices = reshape(X1(1,:,1),3,size(X1,2)/3)';
values = out.LM.permIF;
map = parula(256);
fout = figure();
axes = gca();
colormap(axes,map);
shape.VertexValue = values;
shape.ColorMode = "Indexed";
shape.Material = 'Dull';
shape.ViewMode = 'solid';
shape.RenderAxes = axes;
shape.Visible = true;
%     scan.PatchHandle.FaceColor = 'flat';
axis(axes,'image');
axis(axes,'off');
colorbar(axes,'off');   
% renderBrainSurface(rend,VertexValues{i},fout.ax1{i});
colorbar(axes,'SouthOutside');
shape.VertexSize = 10;
view(axes,1,0);

%%
toCheckSizes= [10, 20, 50, 100];
ret = checkSize(X1, X2, toCheckSizes);
plotExp(ret, toCheckSizes, 'Population Size');

%%
toCheckLSizes= [50, 100, 150, 200, 250];
numExp = length(toCheckLSizes);
permFPsL = zeros(numExp,1);
permDPsL = zeros(numExp,1);
permIPsL = zeros(numExp,1);

for i=1:numExp
    inputSize = toCheckLSizes(i);
    subsample_vec = 3 * (1: ceil((size(X1,2)/(3 * inputSize))):(size(X1,2)/3));
    subsample_vec = ([ (subsample_vec-2)' (subsample_vec -1)' subsample_vec'])';
    subsample_vec = subsample_vec(:);
    ret  = ProcrustesAnova2WayAsymmetryMEM(X1(:, subsample_vec,:),X2(:,subsample_vec,:),200);
    permFPsL(i) = ret.Total.permFF;
    permDPsL(i) = ret.Total.permDF;
    permIPsL(i) = ret.Total.permIF;
end
%%

figure();
plot(toCheckLSizes, permDPsL,'r');
hold on;
plot(toCheckLSizes, permIPsL, 'g');
plot(toCheckLSizes, permFPsL,'b');
ylabel('p-value');
xlabel('landmarks number');
legend("Directional","Individual", "Fluctuating");
title("Dependency of ANOVA2-way asymmetry significance from Number of Landmarks");
hold off;

%%
toCheckRepsNum= [2,4,6];
numExp = length(toCheckRepsNum);
permFPsR = zeros(numExp,1);
permDPsR = zeros(numExp,1);
permIPsR = zeros(numExp,1);

for i=1:numExp
    inputReps = toCheckRepsNum(i);
    subsample_vec = 1: ceil((size(X1,3)/inputReps)):size(X1,3);
    ret  = ProcrustesAnova2WayAsymmetryMEM(X1(:,:,subsample_vec),X2(:,:,subsample_vec),200);
    permFPsR(i) = ret.LM.permFP;
    permDPsR(i) = ret.LM.permDP;
    permIPsR(i) = ret.LM.permIP;
end
%%

figure();
plot(toCheckRepsNum, permDPsR,'r');
hold on;
plot(toCheckRepsNum, permIPsR, 'g');
plot(toCheckRepsNum, permFPsR,'b');
ylabel('p-value');
xlabel('replications number');
legend("Directional","Individual", "Fluctuating");
title("Dependency of ANOVA2-way asymmetry significance from Number of Replications");
hold off;

%%
toCheckPermsNum= [50,100,180,300];
numExp = length(toCheckPermsNum);
permFPsP = zeros(numExp,1);
permDPsP = zeros(numExp,1);
permIPsP = zeros(numExp,1);

for i=1:numExp
    inputReps = toCheckPermsNum(i);
    ret  = ProcrustesAnova2WayAsymmetryMEM(X1,X2,inputReps);
    permFPsP(i) = ret.LM.permFP;
    permDPsP(i) = ret.LM.permDP;
    permIPsP(i) = ret.LM.permIP;
end
%%

figure();
plot(toCheckPermsNum, permFPsP,'r');
hold on;
plot(toCheckPermsNum, permDPsP, 'g');
plot(toCheckPermsNum, permIPsP,'b');
ylabel('p-value');
xlabel('permutations number');
legend("Directional","Individual", "Fluctuating");
title("Dependency of ANOVA2-way asymmetry significance from Number of Permutations");
hold off;

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
f = figure;f.Position = [95  98  2192  1106];f.Color = [1 1 1];%
i=1;
VertexValues{i} = showstruct.LM.I;titlenames{i} = 'I';i=i+1;
val = showstruct.LM.permIF;res = zeros(size(val));
res(val<=0.05) = 0.5;
VertexValues{i} = showstruct.LM.IF;titlenames{i} = 'IF';i=i+1;
res(val<=0.001) = 1;
VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
VertexValues{i} = showstruct.LM.D;titlenames{i} = 'D';i=i+1;
VertexValues{i} = showstruct.LM.DF;titlenames{i} = 'DF';i=i+1;
val = showstruct.LM.permDF;res = zeros(size(val));
res(val<=0.05) = 0.5;%% was 0.002
res(val<=0.001) = 1;
VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
VertexValues{i} = showstruct.LM.F;titlenames{i} = 'F';i=i+1;
VertexValues{i} = showstruct.LM.FF;titlenames{i} = 'FF';i=i+1;
val = showstruct.LM.permFF;res = zeros(size(val));
res(val<=0.05) = 0.5;
res(val<=0.001) = 1;
VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
nValues = length(VertexValues);
arrange = [3 6];
counter = 0;
map = parula(256);
clim = [];
rend = Render{1};
for i=1:nValues
        %i=1;
       counter = counter+1;
       fout.ax1{i} = subplot(arrange(1),arrange(2),counter,'Parent',f);
       colormap(fout.ax1{i},map);
       if ~isempty(clim), set(fout.ax1{i},'clim',clim);end
       renderBrainSurface(rend,VertexValues{i},fout.ax1{i});
       colorbar(fout.ax1{i},'SouthOutside');
       if mod(i,3)>0,          
           set(fout.ax1{i},'clim',[0 max(VertexValues{i})]);
       end
       if mod(i,3)==0,set(fout.ax1{i},'clim',[0 1]); colormap(fout.ax1{i},'summer');end
       view(fout.ax1{i},rend.viewval(1),0);
       light = camlight(fout.ax1{i},'headlight');
       set(light,'Position',get(fout.ax1{i},'CameraPosition'));
       drawnow;
       title(fout.ax1{i},titlenames{i})
       counter = counter+1;
       fout.ax2{i} = subplot(arrange(1),arrange(2),counter,'Parent',f);
       renderBrainSurface(rend,VertexValues{i},fout.ax2{i});
       view(fout.ax2{i},-1*rend.viewval(1),0);
       colorbar(fout.ax2{i},'SouthOutside');
       if mod(i,3)>0,
           set(fout.ax2{i},'clim',[0 max(VertexValues{i})]);
       end
       colormap(fout.ax2{i},map);
       if mod(i,3)==0,set(fout.ax2{i},'clim',[0 1]); colormap(fout.ax2{i},'summer');end
       light = camlight(fout.ax2{i},'headlight');
       set(light,'Position',get(fout.ax2{i},'CameraPosition'));
       drawnow;
       if ~isempty(clim), set(fout.ax2{i},'clim',clim);end
end
% 
% 
% 
% %%
% 
% 
% 
% %save([savepath 'TwoWayProcrustesAnovat0vFinal'],'showstruct','-v7.3');