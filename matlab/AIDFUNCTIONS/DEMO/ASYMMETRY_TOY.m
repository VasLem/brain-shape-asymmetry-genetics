close all;clear;
restoredefaultpath;
% addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));
addpath(genpath('AIDFUNCTIONS'));

file = "../SAMPLE_DATA/augmented_mosquito3d_reps_2.csv";
opts = detectImportOptions(file);
preview(file,opts)
t = readtable(file,opts);

%%
n = size(t,1);
info = zeros(n,3);
for c=1:n
s = string(t.index(c));
p = split(s,'_');
info(c,1) =  str2double(p(1));
info(c,2) = str2double(p(2));
info(c,3) = str2double(p(3));
end
%%
m = readmatrix(file);
m = m(:,2:end);
m = permute(reshape(m,[120, 3, 36]),[3,2,1]);
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
RepShapes = zeros(size(Shapes,1)/2,size(Shapes,2),nRep,'single');
for i=1:nRep
    RepShapes(:,:,i) = Shapes(Sinfo(:,3)==i,:);
end
RepShapesInt16 = int16(RepShapes.*10000);clear RepShapes;
%%
X1 = RepShapesInt16(1:size(RepShapesInt16,1)/2,:,:);
X2 = RepShapesInt16(size(RepShapesInt16,1)/2+1:end,:,:); 
%%
toCheckSizes= [5,8, 10];
ret = checkSize(X1, X2, toCheckSizes,true);  
plotExp(ret, toCheckSizes, 'Population Size');

%%
out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,100);