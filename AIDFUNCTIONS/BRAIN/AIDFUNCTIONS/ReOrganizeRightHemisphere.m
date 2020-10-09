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
IND = 1:length(MASK);
SegmentVal = 1:nSegments;
cmap = colorcube(length(IND));
f = figure;f.Position = [44   592   879   734];
f.Color = [1 1 1];ax = gca;
peer = plotHierarchicalLabelsDGors(ax,rend,SegmentVal,[0 size(cmap,1)], cmap, []);

MASK = rend.UMASK;
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

rend.ReOrd.PerQ.Index = Index;
rend.ReOrd.PerQ.Label = Label;


index = find(Segments.UMASK);
ReOrg.SegIndex = 0*ReOrg.Index;
for i=1:length(index)
    % i=1;
    ReOrg.SegIndex(ReOrg.Index==index(i)) = i;
end
ind = find(ReOrg.SegIndex);
ReOrg.SegIndex = ReOrg.SegIndex(ind);
ReOrg.SegLabel = ReOrg.Label(index);






%% SOME COLOR DEFINITIONS
Color = cell(7,2);
Color{1,1} = [247,247,247]/255;
Color{2,1} = [253,219,199]/255;
Color{3,1} = [209,229,240]/255;
Color{4,1} = [178,24,43]/255;
Color{5,1} = [239,138,98]/255;
Color{6,1} = [103,169,207]/255;
Color{7,1} = [33,102,172]/255;

Color{1,2} = [247,247,247]/255;
Color{2,2} = [209,229,240]/255;
Color{3,2} = [253,219,199]/255;
Color{4,2} = [103,169,207]/255;
Color{5,2} = [33,102,172]/255;
Color{6,2} = [178,24,43]/255;
Color{7,2} = [239,138,98]/255;
%% FACIAL SEGMENTATION
load('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/AM/Shape3DFormat/facerender.mat');%loads facerender
imgpath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/FACEANALYSIS/FIGURES/';
%cmap = colorcube(63);
SegmentVal = 1:63;
HI = facerender.UHI;
cmap = [];
cmap(1,:) = Color{1,1};
cmap(2,:) = Color{2,1};
cmap(3,:) = Color{3,1};
for i=4:1:7
    %i=4;
    cmap(i,:) = Color{i,1};
    ind = sort(getAllChildren(HI,i));
    for j=1:length(ind)
        if mod(j,3)==1
            cmap(ind(j),:) = Color{i,1};
        elseif mod(j,3)==2
            cmap(ind(j),:) = Color{i,1}*0.7;    
        else
            tmp = Color{i,1}*1.3;
            tmp(tmp>1) = 1;
            cmap(ind(j),:) = tmp;
        end
    end
end
f = figure;f.Position = [44   592   879   734];
f.Color = [1 1 1];ax = gca;
peer = plotHierarchicalLabelsDGors(ax,facerender,SegmentVal,[0 63], cmap,[]);
print(f,'-dpng','-r300',[imgpath 'FACIAL SEGMENTS ROSSETTE']);
f = figure;f.Position = [44 887 1143 439];f.Color = [1 1 1];
out = FACEOnSurfaceArrange(f,facerender,SegmentVal,1:6,[2 3],cmap,[1 63]);
print(f,'-dpng','-r300',[imgpath 'FACIAL SEGMENTS SURFACE']);
