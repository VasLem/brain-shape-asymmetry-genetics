close all;clear;
restoredefaultpath;
% addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));
addpath(genpath('AIDFUNCTIONS'));

file = "../SAMPLE_DATA/mosquito3d.csv";
opts = detectImportOptions(file);
preview(file,opts)
t = readtable(file,opts);

%%
n = size(t,1);
info = zeros(n,3);
for c=1:n
s = string(t.index(c));
p = split(s,'_');
info(c,1) =  str2double(p(2));
info(c,2) = str2double(p(3));
info(c,3) = str2double(p(4));
end
%%
m = readmatrix(file);
m = m(:,2:end);
m = permute(reshape(m,[40, 3, 36]),[3,2,1]);
%%

[AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(m,int16.empty(5,0),100,true,'best',true,false);


%% 
LHAligned = AlignedShapes(:,:,info(:,2)==1);
RHAligned = AlignedShapes(:,:,info(:,2)==2);
Shapes = cat(3,LHAligned,RHAligned);

LHInfo = info(info(:,2)==1,:);
RHInfo = info(info(:,2)==2,:);
Sinfo = cat(1,LHInfo,RHInfo);


Shapes = permute(Shapes,[2 1 3]);
Shapes = reshape(Shapes,size(Shapes,1)*size(Shapes,2),size(Shapes,3))';
nRep = 2;
RepShapes = zeros(size(Shapes,1)/2,size(Shapes,2),nRep,'single');% noise injected replications
for i=1:nRep
    RepShapes(:,:,i) = Shapes(Sinfo(:,3)==i,:);
end
RepShapesInt16 = int16(RepShapes.*10000);clear RepShapes;
%%
X1 = RepShapesInt16(1:10,:,:);
X2 = RepShapesInt16(11:end,:,:); 
out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,100);