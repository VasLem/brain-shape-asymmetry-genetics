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
strDATA = cell(1,2);
strRender = cell(1, 2);
savePath = [DATA_DIR 'asym_data.mat'];
disp(['Extracting mat file of serialized objects to ' savePath]);
for r=1:nR
    %r=2
    regphenopath = [phenopath Regions{r} '/'];
    strDATA{r} =  load([regphenopath 'STAGE00DATA']);
    strDATA{r}.Region.AvgShape = strDATA{r}.Region.AvgShape.obj2struc();
    strRender{r} = load([regphenopath 'RENDERMATERIAL.mat']);
    strRender{r}.RefScan = strRender{r}.RefScan.obj2struc();
    
end
save(savePath,  'strDATA', 'strRender');
clear strDATA strRender;

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


%% GPA
% a subselection for now


nSamples = 100;
LH = DATA{1}.Region.AlignedShapes(:, :, 1:nSamples);
RH = DATA{2}.Region.AlignedShapes(:, :, 1:nSamples);
RH(:,1,:,:) = -1*RH(:,1,:,:);
TotalShapes = cat(3,LH,RH);


Template = clone(DATA{1}.Region.AvgShape);
strTemplate = Template.obj2struc();
savePath = [DATA_DIR 'gpa_data.mat'];
save(savePath,  'TotalShapes', 'strTemplate');
[AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(TotalShapes,Template,3,true,false,true,false);
savePath = [DATA_DIR 'aligned.mat'];
save(savePath,  'AlignedShapes');
%%
clear DATA LH RH TotalShapes;
%% TWO WAY PROCRUSTES ANOVA ON REDUCED DATA
% nAnovaSamples = size(AlignedShapes,3)/2;
LHAligned = AlignedShapes(:,:,1:nSamples);% a subselection for now
RHAligned = AlignedShapes(:,:,nSamples+1:nSamples+nSamples);
% RHAligned = reshape(permute(repmat(1:(size(RHAligned, 1) * size(RHAligned, 2)), nSamples, 1), [2, 1]), size(RHAligned))/prod(size(RHAligned, [1,2]));
% LHAligned = reshape(permute(-repmat(1:(size(LHAligned, 1) * size(LHAligned, 2)), nSamples, 1), [2, 1]), size(LHAligned))/prod(size(LHAligned, [1,2]));
    %%
Shapes = cat(3,LHAligned,RHAligned);
% Shapes = Shapes(1:100:end,:,:);% reducing the amount of vertices
Shapes = permute(Shapes,[2 1 3]);
Shapes = reshape(Shapes,size(Shapes,1)*size(Shapes,2),size(Shapes,3))';
%%
% nRep = 3
% RepShapes = zeros(size(Shapes,1),size(Shapes,2),nRep,'single');% noise injected replications
% for i=1:nRep
%     RepShapes(:,:,i) = single(Shapes) + single(randn(size(Shapes,1),size(Shapes,2)).*0.05);
% end
nRep = 3;
RepShapes = zeros(size(Shapes,1),size(Shapes,2),nRep,'single');% noise injected replications
for i=1:nRep
    % RepShapes(:,:,i) = single(Shapes);
    RepShapes(:,:,i) = single(Shapes) + single(randn(size(Shapes,1),size(Shapes,2)).*0.05);
end
%%
RepShapesInt16 = int16(RepShapes.*10000);clear RepShapes;
%%
X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:);
totalX = cat(1, X1, X2);
savePath = [DATA_DIR 'replicated.mat'];
save(savePath,  'totalX');
%%
out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,nSamples);
%% BELOW IS AN IDEA OF RENDERING, BUT WILL NOT WORK BECAUSE WE DO NOT HAVE ALL THE MESH POINTS
f = figure;f.Position = [95  98  2192  1106];f.Color = [1 1 1];%
i=1;
VertexValues{i} = out.LM.I;titlenames{i} = 'I';i=i+1;
VertexValues{i} = out.LM.IF;titlenames{i} = 'IF';i=i+1;
val = out.LM.permIF;res = zeros(size(val));
res(val<=0.05) = 0.5;
res(val<=0.001) = 1;
VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
VertexValues{i} = out.LM.D;titlenames{i} = 'D';i=i+1;
VertexValues{i} = out.LM.DF;titlenames{i} = 'DF';i=i+1;
val = out.LM.permDF;res = zeros(size(val));
res(val<=0.002) = 0.5;
res(val<=0.001) = 1;
VertexValues{i} = res;titlenames{i} = 'p';i=i+1;
VertexValues{i} = out.LM.F;titlenames{i} = 'F';i=i+1;
VertexValues{i} = out.LM.FF;titlenames{i} = 'FF';i=i+1;
val = out.LM.permFF;res = zeros(size(val));
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



%%



%save([savepath 'TwoWayProcrustesAnovat0vFinal'],'out','-v7.3');
%% TWO WAY PROCRUSTES ANOVA
nSamples = size(AlignedShapes,3)/2;
LHAligned = AlignedShapes(:,:,1:nSamples);
RHAligned = AlignedShapes(:,:,nSamples+1);
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
save([savepath 'TwoWayProcrustesAnovat0vFinal'],'out','-v7.3');
%% INVESTIGATE NUMBER OF SAMPLES ON DA

in = load([savepath 'TwoWayProcrustesAnovat0v4']);
in5 = load([savepath 'TwoWayProcrustesAnovat0v5']);

n = 100:100:500;
N = length(n);
TOT = cell(1,N);
avgDA = zeros(1,N);

for i=1:N
    disp(num2str(i));
    X1 = RepShapesInt16(1:n(i),:,:);
    X2 = RepShapesInt16(nSamples+1:nSamples+n(i),:,:);
    out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,0);
    TOT{i} = out;
end
save([savepath 'SAMPLESIZEINFL'],'TOT','-v7.3');


avgDA = zeros(1,N);
for i=1:N
    TOT{i}.Total
    avgDA(i) = TOT{i}.Total.D;
end
    




%%

addpath('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/CODE/ANALYSIS/');
addpath('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/CODE/ANALYSIS/AIDFUNCTIONS/');
load([savepath 'TwoWayProcrustesAnovat0vFinal']);
f = figure;f.Position = [95         261        1185         943];f.Color = [1 1 1];%
i=1;
VertexValues = [];
VertexValues{i} = out.LM.I;titlenames{i} = 'I';i=i+1;
VertexValues{i} = out.LM.IF;titlenames{i} = 'IF';i=i+1;
%VertexValues{i} = out.LM.IP<=0.05;titlenames{i} = 'p';i=i+1;
VertexValues{i} = out.LM.D;titlenames{i} = 'D';i=i+1;
VertexValues{i} = out.LM.DF;titlenames{i} = 'DF';i=i+1;
%VertexValues{i} = out.LM.DP<=0.05;titlenames{i} = 'p';i=i+1;
VertexValues{i} = out.LM.F;titlenames{i} = 'F';i=i+1;
VertexValues{i} = out.LM.FF;titlenames{i} = 'FF';i=i+1;
%VertexValues{i} = out.LM.FP<=1e-5;titlenames{i} = 'p';i=i+1;
nValues = length(VertexValues);
arrange = [3 4];
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
       %if mod(i,2)>0
           colorbar(fout.ax1{i},'SouthOutside');
           set(fout.ax1{i},'clim',[0 max(VertexValues{i})]);
       %end
       %if mod(i,3)==0,set(fout.ax1{i},'clim',[0 1]);end
       if mod(i,2) == 1, set(fout.ax1{i},'clim',[0 0.03]); end
       view(fout.ax1{i},rend.viewval(1),0);
       light = camlight(fout.ax1{i},'headlight');
       set(light,'Position',get(fout.ax1{i},'CameraPosition'));
       drawnow;
       %title(fout.ax1{i},titlenames{i})
       counter = counter+1;
       fout.ax2{i} = subplot(arrange(1),arrange(2),counter,'Parent',f);
       renderBrainSurface(rend,VertexValues{i},fout.ax2{i});
       view(fout.ax2{i},-1*rend.viewval(1),0);
       %if mod(i,3)>0,
           colorbar(fout.ax2{i},'SouthOutside');
           set(fout.ax2{i},'clim',[0 max(VertexValues{i})]);
       %end
       %if mod(i,3)==0,set(fout.ax2{i},'clim',[0 1]);end
       if mod(i,2) == 1, set(fout.ax2{i},'clim',[0 0.03]); end
       light = camlight(fout.ax2{i},'headlight');
       set(light,'Position',get(fout.ax2{i},'CameraPosition'));
       drawnow;
       colormap(fout.ax2{i},map);
       if ~isempty(clim), set(fout.ax2{i},'clim',clim);end
end
print(f,'-dpng','-r300',[imgpath 'TwoWayProcrustesAnovaFinal']);

%%

% imaging outcome

imgpath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/ASYMMETRY/';


left = Render{1}.RefScan;
right = Render{2}.RefScan;

v = viewer(left)
viewer(right,v);

v2 = viewer(right);


reflright = clone(right);
reflright.Vertices(:,1) = reflright.Vertices(:,1)*-1;

viewer(reflright,v);

Floating = pairWiseSuperimposition(left,reflright);

v1 = viewer(left);
viewer(Floating,v1);

dist = sqrt(sum((left.Vertices-Floating.Vertices).^2,2));

left.VertexValue = dist;
v2 = viewer(left);

f=figure;f.Color = [1 1 1];f.Position = [4 780  1034  546];
ax{1} = subplot(1,2,1);ax{2} = subplot(1,2,2);
viewvals = [-90 90];
cmapH = parula(256);
for a=1:2
    scan = clone(Render{1}.RefScan);
    scan.VertexValue = dist;
    scan.ColorMode = "Indexed";
    scan.Material = 'Dull';
    scan.ViewMode = 'solid';
    scan.RenderAxes = ax{a};
    scan.Visible = true;
    view(ax{a},viewvals(a),0);
    axis(ax{a},'image');axis(ax{a},'off');colorbar(ax{a},'SouthOutside');   
    light = camlight(ax{a},'headlight');set(light,'Position',get(ax{a},'CameraPosition'));
    colormap(ax{a},cmapH);
end
print(f,'-dpng','-r300',[imgpath 'LeftRightAsymmetryMap']);

%% TO REFINE
