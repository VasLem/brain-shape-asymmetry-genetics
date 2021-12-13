clc;clear all;
addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));
load('/uz/data/avalok/mic/tmp/myuan0/BrainThickness/CODE/Pheno/DATA/COVThickness.mat'); %resT
%outpath = '/uz/data/avalok/mic/tmp/myuan0/BrainThickness/CODE/Pheno/PHENOTYPES60/';

load('/uz/data/avalok/mic/tmp/myuan0/BrainThickness/CODE/Pheno/Segmentation/WeissSegmentation_TS_SPerc80_nL9.mat');
outpath = '/uz/data/avalok/mic/tmp/myuan0/BrainThickness/CODE/Pheno/PHENOTYPES/PHENOTYPES_TS_SPerc80_PercVar60/';
LABELS = LABELS_TS_SPerc80;
MASK = MASK_TS_SPerc80;

%%
index = 1:9; 
u_levels = index(end);% 9 levels to analyze
ULABELS = LABELS(index,:);
UHI = HierarchicalInterface;
UHI.nL = u_levels;
UMASK = MASK(1:UHI.nLC);
disp(['Number of total clusters: ' num2str(length(find(UMASK)))]);
[nL,nV] = size(ULABELS);
nS = size(resT,1);
THICKNESS = resT';
for i=1:UHI.nLC
    disp([num2str(i) '/' num2str(UHI.nLC)]);
    if UMASK(i)~=0
        [l,c] = Ind2LC(UHI,i);
        disp(['Perform PCA for Level ' num2str(l) ' Cluster ' num2str(c)]);
        CLInd = ULABELS(l,:);
        SubInd = find(CLInd==c);
        subTHICKNESS = THICKNESS(SubInd,:);
        tPCA = {};
        tPCA.LC = [l,c];
        tPCA.AvgVec = mean(subTHICKNESS,2);
        D = (double(subTHICKNESS) - repmat(tPCA.AvgVec,1,size(subTHICKNESS,2)))'; % mean-centering
        [tPCA.EigVec,tPCA.Tcoeff,tPCA.EigVal] = my_princomp(D);
        %tPCA.EigStd = sqrt(tPCA.EigVal);
        tPCA.Explained = 100*tPCA.EigVal/sum(tPCA.EigVal);
        percvar = 60; % retain PCs for 60% total variance % <500 PCs for 1st segment
        nr = 1;perc = tPCA.Explained(1);
        while perc < percvar
            perc = perc + tPCA.Explained(nr+1);
            nr = nr+1;
        end 
        tPCA.percPA = percvar; 
        tPCA.nrPA = nr;                       % Number of Principal Components
        tPCA.EigVal = tPCA.EigVal(1:nr);      % Eigen Values
        tPCA.EigStd = sqrt(tPCA.EigVal);      % Standard deviations per component
        tPCA.EigVec = tPCA.EigVec(:,1:nr);    % Eigen Vectors
        tPCA.Tcoeff = tPCA.Tcoeff(:,1:nr);    % myPCA coeff of training data
        tPCA.Explained = tPCA.Explained(1:nr);% Percentage Explained in Eigenvectors
        tPCA.n = size(tPCA.Tcoeff,1);         % number of trainings Data
        tPCA.DepVar = tPCA.Tcoeff./repmat(tPCA.EigStd',tPCA.n,1);
        save([outpath 'tPC_Ind' num2str(i) '_L' num2str(l) '_C' num2str(c)],'tPCA','-v7.3');    
    end
end
