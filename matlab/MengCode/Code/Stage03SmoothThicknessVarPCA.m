clc;clear all;
addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));

path = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainThickness/CODE/Pheno';
load([path '/DATA/COVThickness.mat']); %resT
load([path '/Segmentation/WeissSegmentation_Tsmooth_iter30_nL9.mat']);
outpath = [path '/PHENOTYPES/PHENOTYPES_smoothT_iter30_PercVar60/'];

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
            %SPACE = shapePCA;
            %scan = crop(RefScan,'VertexIndex',SubInd);
            %SPACE.RefScan = clone(scan);
            %getAverage(SPACE,redTHICKNESS); 
            %getModel(SPACE,redTHICKNESS); 
            %D = LSGenProcrustes(SPACE,redTHICKNESS,true,3,scan); 
            %getAverage(SPACE,D);
            %getModel(SPACE,D); 
            % ALTERNATIVE is fixed and low percentage of variance
            % explained: 
            %percPA = 80;
            %stripPercVar(SPACE,percPA);
            %DepVar = SPACE.Tcoeff./repmat(SPACE.EigStd',SPACE.n,1);
           % nrPA = SPACE.nrEV;
            %SPACE = clone(SPACE);
         save([outpath 'tPC_Ind' num2str(i) '_L' num2str(l) '_C' num2str(c)],'tPCA','-v7.3');    
    end
end

% % check compactness
% outpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainThickness/CODE/Pheno/PHENOTYPES/PHENOTYPES_T_Dist80_PercVar60/';
% % outpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainThickness/CODE/Pheno/PHENOTYPES/PHENOTYPES_TS_SPerc20_PercVar60/';
% 
% cd(outpath);files = dir('**/*.mat');
% nfiles = length(files);
% compact_score = zeros(nfiles,3);
% for i= 1:nfiles
%     in = load([outpath files(i).name]);  
%     fpart = split(files(i).name,["_Ind","_L"]); 
%     compact_score(i,1) = str2double(fpart{2,1}); % first column for Ind
%     compact_score(i,2) = in.tPCA.nrPA;
%     compact_score(i,3) = in.tPCA.nrPA/in.tPCA.percPA;
% end
% compact_score = sortrows(compact_score);
% mean(compact_score(1:50,:))
% mean(compact_score)

