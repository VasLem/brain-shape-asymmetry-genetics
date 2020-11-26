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
Ns = 100;
LH = DATA{1}.Region.AlignedShapes(1:Ns:end, :, 1:nSamples);
RH = DATA{2}.Region.AlignedShapes(1:Ns:end, :, 1:nSamples);
RH(:,1,:,:) = -1*RH(:,1,:,:);
Template = clone(DATA{1}.Region.AvgShape);
reducedTemplate = shape3D;
reducedTemplate.struc2obj(Template.obj2struc());

reducedTemplate.Vertices = reducedTemplate.Vertices(1:Ns:end,:);
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



nRep = 3;
% 
input_to_r_path =  [DATA_DIR 'm2r.mat'];
save(input_to_r_path, 'LH', 'RH','strReducedTemplate', 'strTemplate', 'nSamples','Ns','nRep');

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

display3DLandmarks(reducedTemplate, reducedTemplateAligned_R);
%%
display3DLandmarks(LHAligned(:,:,1), LHAligned_R(:,:,1));

%%
display3DLandmarks(reducedTemplateAligned_R, LHAligned_R(:,:,1));
%%

display3DLandmarks(LHAligned(:,:,1), reducedTemplate.Vertices);
%%


% figure;
% hold on;
% scatter3(RHAligned(:,1,1),RHAligned(:,2,1),RHAligned(:,3,1),'b');,
% scatter3(LHAligned (:,1,1),LHAligned (:,2,1),LHAligned (:,3,1),'r');,
% hold off;


LHAlignedInt16 = int16(LHAligned.*10000);clear LHAligned;
RHAlignedInt16 = int16(RHAligned.*10000);clear RHAligned;
%%

Shapes = permute(AlignedShapes,[2 1 3]);
Shapes = reshape(Shapes,size(Shapes,1)*size(Shapes,2),size(Shapes,3))';
clear AlignedShapes;
%nRep = 6;
nRep = 3;
RepShapes = zeros(size(Shapes,1),size(Shapes,2),nRep,'single');
for i=1:1:nRep
    RepShapes(:,:,i) = single(Shapes) + single(randn(size(Shapes,1),size(Shapes,2)).*0.05);
end
RepShapesInt16 = int16(RepShapes.*10000);clear RepShapes;

X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:);

out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,0);




%% BELOW IS AN IDEA OF RENDERING, BUT WILL NOT WORK BECAUSE WE DO NOT HAVE ALL THE MESH POINTS
% f = figure;f.Position = [95  98  2192  1106];f.Color = [1 1 1];%
% i=1;
% VertexValues{i} = out.LM.I;titlenames{i} = 'I';i=i+1;
% val = out.LM.permIF;res = zeros(size(val));
% res(val<=0.05) = 0.5;
% VertexValues{i} = out.LM.IF;titlenames{i} = 'IF';i=i+1;
% res(val<=0.001) = 1;
% VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
% VertexValues{i} = out.LM.D;titlenames{i} = 'D';i=i+1;
% VertexValues{i} = out.LM.DF;titlenames{i} = 'DF';i=i+1;
% val = out.LM.permDF;res = zeros(size(val));
% res(val<=0.002) = 0.5;
% res(val<=0.001) = 1;
% VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
% VertexValues{i} = out.LM.F;titlenames{i} = 'F';i=i+1;
% VertexValues{i} = out.LM.FF;titlenames{i} = 'FF';i=i+1;
% val = out.LM.permFF;res = zeros(size(val));
% res(val<=0.05) = 0.5;
% res(val<=0.001) = 1;
% VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
% nValues = length(VertexValues);
% arrange = [3 6];
% counter = 0;
% map = parula(256);
% clim = [];
% rend = Render{1};
% for i=1:nValues
%         %i=1;
%        counter = counter+1;
%        fout.ax1{i} = subplot(arrange(1),arrange(2),counter,'Parent',f);
%        colormap(fout.ax1{i},map);
%        if ~isempty(clim), set(fout.ax1{i},'clim',clim);end
%        renderBrainSurface(rend,VertexValues{i},fout.ax1{i});
%        colorbar(fout.ax1{i},'SouthOutside');
%        if mod(i,3)>0,          
%            set(fout.ax1{i},'clim',[0 max(VertexValues{i})]);
%        end
%        if mod(i,3)==0,set(fout.ax1{i},'clim',[0 1]); colormap(fout.ax1{i},'summer');end
%        view(fout.ax1{i},rend.viewval(1),0);
%        light = camlight(fout.ax1{i},'headlight');
%        set(light,'Position',get(fout.ax1{i},'CameraPosition'));
%        drawnow;
%        title(fout.ax1{i},titlenames{i})
%        counter = counter+1;
%        fout.ax2{i} = subplot(arrange(1),arrange(2),counter,'Parent',f);
%        renderBrainSurface(rend,VertexValues{i},fout.ax2{i});
%        view(fout.ax2{i},-1*rend.viewval(1),0);
%        colorbar(fout.ax2{i},'SouthOutside');
%        if mod(i,3)>0,
%            set(fout.ax2{i},'clim',[0 max(VertexValues{i})]);
%        end
%        colormap(fout.ax2{i},map);
%        if mod(i,3)==0,set(fout.ax2{i},'clim',[0 1]); colormap(fout.ax2{i},'summer');end
%        light = camlight(fout.ax2{i},'headlight');
%        set(light,'Position',get(fout.ax2{i},'CameraPosition'));
%        drawnow;
%        if ~isempty(clim), set(fout.ax2{i},'clim',clim);end
% end
% 
% 
% 
% %%
% 
% 
% 
% %save([savepath 'TwoWayProcrustesAnovat0vFinal'],'out','-v7.3');