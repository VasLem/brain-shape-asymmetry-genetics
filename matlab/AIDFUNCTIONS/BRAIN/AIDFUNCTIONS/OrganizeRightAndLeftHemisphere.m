%% LOADING BRAIN GWAS
%close all; clear all;
addpath('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/CODE/ANALYSIS/');
addpath('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/CODE/ANALYSIS/AIDFUNCTIONS');
studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/';
cd(studypath);
%% PHENO INFO
in = load('/uz/data/avalok/mic/tmp/SHARED/pclaes4/BRAIN/UKBIOBANK/PHENOTYPES/LH/AUX_LH_PA80_20190928.mat');
IID = in.DB.IID;clear in;
%% LOADING RENDER MATERIAL
Render{1} = load('/home/pclaes4/AVALOKTMP/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/LH/RENDERMATERIAL.mat');
Render{2} = load('/home/pclaes4/AVALOKTMP/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/RH/RENDERMATERIAL.mat');
in  = load('/home/pclaes4/AVALOKTMP/==MATLAB==/myToolboxes/COLORMAPS/EffectSizeColorMap2.mat');
cmapP = in.EffectSizeColorMap2;
%% TRY AND REORDER RIGHT HEMISPHERE
rend = Render{2};
nSegments = length(find(rend.UMASK));
IND = 1:length(rend.UMASK);
SegmentVal = 1:nSegments;
cmap = colorcube(length(IND));
f = figure;f.Position = [44   592   879   734];
f.Color = [1 1 1];ax = gca;
peer = plotHierarchicalLabelsDGors(ax,rend,SegmentVal,[0 size(cmap,1)], cmap, []);

HI = rend.UHI;
ind = [];
shift = [0 1 2 4 8 16 32 64 128];
for l=1:HI.nL
    tmp = [];
    cl = 1:2^(l-1);
    for c=1:length(cl)
        tmp = [tmp LC2Ind(HI,l,c)];
    end
    tmp = circshift(tmp,shift(l));
    ind = [ind tmp];
end
rend.ReOrder.Ind = ind;
rend.ReOrder.UMASK = rend.UMASK(rend.ReOrder.Ind);
rend.ReOrder.UHI = rend.UHI;

f = figure;f.Position = [44   592   879   734];
f.Color = [1 1 1];ax = gca;
out = reOrderRightSegmentValues(SegmentVal,rend); 
peer = plotHierarchicalLabelsDGors(ax,rend.ReOrder,out,[0 size(cmap,1)], cmap, []);
LABELS = rend.ULABELS;
for i=1:HI.nLC
   [l,c] = Ind2LC(HI,i);
   [oldl,oldc] = Ind2LC(HI,rend.ReOrder.Ind(i));
   if ~(l==oldl), error('levels do not match');end
   LABELS(l,find(rend.ULABELS(l,:)==oldc)) = c;
end
rend.ReOrder.ULABELS = LABELS;
rend.ReOrder.viewval = 90;
rend.ReOrder.RefScan = clone(rend.RefScan);


f = figure;f.Position = [44 887 1143 439];f.Color = [1 1 1];
BRAINOnSurfaceArrangev2(f,rend,SegmentVal,1:9,[3 3],cmap,[1 size(cmap,1)]);

f = figure;f.Position = [44 887 1143 439];f.Color = [1 1 1];
BRAINOnSurfaceArrangev2(f,rend.ReOrder,out,1:9,[3 3],cmap,[1 size(cmap,1)]);
% NOW ADDING THE REORGANIZATION INDEXING PER QUADRANT

Index = [1 2 3 4 sort(getAllChildren(rend.UHI,4)) 5 sort(getAllChildren(rend.UHI,5)) 6 sort(getAllChildren(rend.UHI,6)) 7 sort(getAllChildren(rend.UHI,7))];
Label = [1 2 3 4*ones(1,length([4 sort(getAllChildren(rend.UHI,4))])) 5*ones(1,length([5 sort(getAllChildren(rend.UHI,5))])) 6*ones(1,length([6 sort(getAllChildren(rend.UHI,6))])) 7*ones(1,length([7 sort(getAllChildren(rend.UHI,7))]))];

rend.PerQ.Index = Index;
rend.PerQ.Label = Label;

rend.ReOrder.PerQ.Index = Index;
rend.ReOrder.PerQ.Label = Label;

index = find(rend.UMASK);
SegIndex = 0*Index;
for i=1:length(index)
    % i=1;
    SegIndex(Index==index(i)) = i;
end
ind = find(SegIndex);
rend.PerQ.SegIndex = SegIndex(ind);
rend.PerQ.SegLabel = Label(index);

index = find(rend.ReOrder.UMASK);
SegIndex = 0*Index;
for i=1:length(index)
    % i=1;
    SegIndex(Index==index(i)) = i;
end
ind = find(SegIndex);
rend.ReOrder.PerQ.SegIndex = SegIndex(ind);
rend.ReOrder.PerQ.SegLabel = Label(index);
save('/home/pclaes4/AVALOKTMP/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/RH/RENDERMATERIALREORDERED.mat','rend');
%% REORDER LEFT HEMISHPERE
rend = Render{1};
% nothing changes for the left hemisphere
rend.ReOrder.Ind = 1:length(rend.UMASK);
rend.ReOrder.UMASK = rend.UMASK;
rend.ReOrder.UHI = rend.UHI;
rend.ReOrder.ULABELS = rend.ULABELS;
rend.ReOrder.viewval = rend.viewval;
rend.ReOrder.RefScan = clone(rend.RefScan);

% NOW ORGANIZING PER QUADRANT

Index = [1 2 3 4 sort(getAllChildren(rend.UHI,4)) 5 sort(getAllChildren(rend.UHI,5)) 6 sort(getAllChildren(rend.UHI,6)) 7 sort(getAllChildren(rend.UHI,7))];
Label = [1 2 3 4*ones(1,length([4 sort(getAllChildren(rend.UHI,4))])) 5*ones(1,length([5 sort(getAllChildren(rend.UHI,5))])) 6*ones(1,length([6 sort(getAllChildren(rend.UHI,6))])) 7*ones(1,length([7 sort(getAllChildren(rend.UHI,7))]))];

rend.PerQ.Index = Index;
rend.PerQ.Label = Label;

rend.ReOrder.PerQ.Index = Index;
rend.ReOrder.PerQ.Label = Label;

index = find(rend.UMASK);
SegIndex = 0*Index;
for i=1:length(index)
    % i=1;
    SegIndex(Index==index(i)) = i;
end
ind = find(SegIndex);
rend.PerQ.SegIndex = SegIndex(ind);
rend.PerQ.SegLabel = Label(index);

index = find(rend.ReOrder.UMASK);
SegIndex = 0*Index;
for i=1:length(index)
    % i=1;
    SegIndex(Index==index(i)) = i;
end
ind = find(SegIndex);
rend.ReOrder.PerQ.SegIndex = SegIndex(ind);
rend.ReOrder.PerQ.SegLabel = Label(index);

save('/home/pclaes4/AVALOKTMP/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/LH/RENDERMATERIALREORDERED.mat','rend');
%% THE END
