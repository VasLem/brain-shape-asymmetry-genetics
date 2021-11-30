%%
close all;clear all;
restoredefaultpath;
addpath(genpath('AIDFUNCTIONS'));

DATA_DIR = '../SAMPLE_DATA/';
PHENO_DIR = [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/PHENOTYPES/'];
primarycovpath = '/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/COVARIATES';
secondarycovpath = '../SAMPLE_DATA/IMAGEN/BRAIN/BATCH2_DATA/';
% covpath = [datapath 'COVARIATES/'];
% phenopath = [datapath 'PHENOTYPING/'];
% load([phenopath 'FinalImageIDs_20407subjs']);%loads imageIDs
% IMGID = num2cell(imageIDs);
% IMGID = listofnum2str(IMGID);
% cd(covpath);

%% GETTING RELEVANT DATA ENTRY FIELDS
x = readtable(fullfile(primarycovpath,'Final_list_covariates_GWAS.xlsx' ));
COVIDs = x.Var1;
CovNames = x.Var2;
nCOV = length(COVIDs);
%%

ids = load(fullfile(PHENO_DIR, 'BATCH2_2021_DATA_IID')).iids;
csvpath = fullfile(primarycovpath, 'ukb35648.csv');
[parsedData1, ids1, cols1] = readcsv(csvpath, ids, COVIDs);

csvpath = fullfile(secondarycovpath, 'ukb44263', 'ukb44263.csv');
[parsedData2, ids2, cols2] = readcsv(csvpath, ids, COVIDs);
%% TODO, adapt the following to connect the two datasets, check ids first!
for c = 1:length(uniqueInds)
    subset = data(:, colInds == c);
    ret = zeros(length(mapped_to), 1);
    ret(:) = nan;
    ret = num2cell(ret);
    for col = 1:size(subset,2)
        mask = ~cellfun(@isempty, subset(:,col));
        ret(mask) = subset(mask, col);
    end
    parsedData(:, c) = ret;
end


function [parsedData, filteredIds, filteredCols ]= readcsv(csvpath, ids, COVIDs)
fid = fopen(csvpath);
headers = strsplit(fgetl(fid), ',');
fclose(fid);
colIndices = arrayfun(@(x) find(startsWith(headers, ['"', num2str(x) '-'])), COVIDs , 'UniformOutput', false);
found = ~cellfun(@isempty, colIndices);
fCOVIDs = COVIDs(found);

filteredCols = fCOVIDs;
fCOVIDsInds = colIndices(found);
%%
%%
%%
numCol = length(headers);
fmt_parts = repmat({'%*s'}, 1, numCol);
fmt_parts(horzcat(fCOVIDsInds{:})) = {'%q'};
fmt_parts(1) = {'%q'};
fmt = strjoin(fmt_parts, '');
fid = fopen(csvpath);
data = textscan(fid, fmt, 'Delimiter', ',', 'Headerlines', 1);
data{1} = num2cell(str2double(data{1}));
fclose(fid);
%%
num = 1:length(fCOVIDsInds);
colInds = cellfun(@(x, ind)(repmat(ind, length(x), 1)), fCOVIDsInds, num2cell(num + 1)', 'UniformOutput', false);
colInds = [1; vertcat(colInds{:})];
%%
data = [data{:}];

%%
all_ids =[data{:, 1}];
assignmentMatrix = ( str2double(ids) == all_ids);
[to_keep, mapped_to] = find(assignmentMatrix);
filteredIds = ids(to_keep);
%%
clear assignmentMatrix
%%
data = data(mapped_to, :);
%%
uniqueInds = unique(colInds);
parsedData = zeros(length(mapped_to), length(uniqueInds));
parsedData(:) = nan;
parsedData = num2cell(parsedData);
for c = 1:length(uniqueInds)
    subset = data(:, colInds == c);
    ret = zeros(length(mapped_to), 1);
    ret(:) = nan;
    ret = num2cell(ret);
    for col = 1:size(subset,2)
        mask = ~cellfun(@isempty, subset(:,col));
        ret(mask) = subset(mask, col);
    end
    parsedData(:, c) = ret;
end
end
%%



%%
%
%

% %% FIRST SMALLER CSV FILE
% inD = importdata(fullfile(primarycovpath, 'ukb25183.csv'));
% D = inD;
% headers = D{1};
% step 1A: match image ids 500k->20k subjects
% nD = length(D);
% index = zeros(1,nD);
% for i=1:1:nD
%     if mod(i,10000)==0, disp(num2str(i));end
%     str = D{i};
%     ind = strfind(str,'"');
%     id = str(ind(1)+1:ind(2)-1);
%     ind12 = vlookupFast(IMGID,id);
%     if ~isempty(ind12),index(i) = 1; end
% end
% % % step 1B: reduce the DATA
% keep = find(index);
% D = D(keep); %#ok<*FNDSB>
% nIMG = length(D);
% step 2A: match the COV IDs
%%
% indquotes = strfind(headers,'"');
% nHeaders = length(indquotes)/2;
% index = zeros(2,nCOV);
% indexCOV = zeros(1,nCOV);
% for i=1:1:nCOV
%     %i=1;
%     %ind = strfind(headers,['"' COVIDs{i} '-0.']);
%     ind = strfind(headers,['"' num2str(COVIDs(i)) '-']);
%     if isempty(ind), continue; end
%     indexCOV(i) = 1;
%     listindex = find(indquotes==ind(1));
%     index(1,i) = listindex;
%     index(2,i) = listindex+1;
% end
% % step 2B: reduce the data and store
% keep = find(indexCOV);
% index = index(:,keep);
% COVDATA.nCOV = length(keep);
% COVDATA.COVIDs = COVIDs(keep);
% COVDATA.COVNames = COVNames(keep);
% COVDATA.nSUB = nIMG;
% COVDATA.IMGID = cell(1,nIMG);
% COVDATA.DATA = cell(nIMG,COVDATA.nCOV);
% for i=1:1:nIMG
%     if mod(i,100)==0, disp(num2str(i));end
%     str = D{i};
%     ind = strfind(str,'"');
%     id = str(ind(1)+1:ind(2)-1);
%     COVDATA.IMGID{i} = id;
%     for j=1:1:COVDATA.nCOV
%         COVDATA.DATA{i,j} = str(ind(index(1,j))+1:ind(index(2,j))-1);
%     end
% end
% COVDATA1 = COVDATA;
% save('COVDATA1','COVDATA1','-v7.3');
% %% SECOND CSV FILE
% inD = importdata('ukb35648.csv');
% D = inD;
% headers = D{1};
% % step 1A: match image ids 500k->20k subjects
% nD = length(D);
% index = zeros(1,nD);
% for i=1:1:nD
%     if mod(i,10000)==0, disp(num2str(i));end
%     str = D{i};
%     ind = strfind(str,'"');
%     id = str(ind(1)+1:ind(2)-1);
%     ind12 = vlookupFast(IMGID,id);
%     if ~isempty(ind12),index(i) = 1; end
% end
% % step 1B: reduce the DATA
% keep = find(index);
% D = D(keep); %#ok<*FNDSB>
% nIMG = length(D);
% % step 2A: match the COV IDs
% indquotes = strfind(headers,'"');
% nHeaders = length(indquotes)/2;
% index = zeros(2,nCOV);
% indexCOV = zeros(1,nCOV);
% for i=1:1:nCOV
%     %i=1;
%     ind = strfind(headers,['"' COVIDs{i} '-']);
%     %ind = strfind(headers,['"' '23226' '-']);
%     if isempty(ind), continue; end
%     indexCOV(i) = 1;
%     listindex = find(indquotes==ind(1));
%     index(1,i) = listindex;
%     index(2,i) = listindex+1;
% end
% % step 2B: reduce the data and store
% keep = find(indexCOV);
% index = index(:,keep);
% COVDATA.nCOV = length(keep);
% COVDATA.COVIDs = COVIDs(keep);
% COVDATA.COVNames = COVNames(keep);
% COVDATA.nSUB = nIMG;
% COVDATA.IMGID = cell(1,nIMG);
% COVDATA.DATA = cell(nIMG,COVDATA.nCOV);
% for i=1:1:nIMG
%     if mod(i,100)==0, disp(num2str(i));end
%     str = D{i};
%     ind = strfind(str,'"');
%     id = str(ind(1)+1:ind(2)-1);
%     COVDATA.IMGID{i} = id;
%     for j=1:1:COVDATA.nCOV
%         COVDATA.DATA{i,j} = str(ind(index(1,j))+1:ind(index(2,j))-1);
%     end
% end
% COVDATA2 = COVDATA;
% save('COVDATA2','COVDATA2','-v7.3');
% %% NOTE TO DISCUSS WITH ALEJANDRA
% % Some variables have a lot of missing values, especially the headbone
% % density measures, In the nature paper however, it seems that they do not
% % put these in the final covariate model. I'm not convinced that their
% % reasoning is ok, anything that is associated to brain shape can also be
% % associated to the same genetics, that is the nature of the game, it does
% % not mean your association is false, it simply means your genetic
% % functionality is not known. It is only a real problem, if these variables
% % associated with almost all the genome (like ancestry and stratification
% % does). But from Supplementary figure 22, this does not seem to be the
% % case
% %% IMPORTING GENETIC AXIS
% in = load('UKB_EUR_19670_PCs');
% nPC = size(in.AA,2);
% COVDATA3.IID = in.IID;
% COVDATA3.FID = in.FID;
% COVDATA3.COVIDs = cell(nPC,1);
% COVDATA3.COVNames = cell(nPC,1);
% for i=1:nPC
%     COVDATA3.COVIDs{i} = ['00' num2str(i)];
%     COVDATA3.COVNames{i} = ['KUL GENO PC ' num2str(i)];
% end
% COVDATA3.DATA = num2cell(in.AA);
% save('COVDATA3','COVDATA3','-v7.3');
% %% ASSEMBLING EVERYTHING
% COVDATA = COVDATA3;
% COVDATA.nSUB = length(COVDATA.IID);
% % adding first set of covariates
% addCOVDATA = COVDATA1;
% [ind12,~] = vlookupFast(COVDATA.IID,addCOVDATA.IMGID);
% index = find(~ismember(addCOVDATA.COVIDs,COVDATA.COVIDs));
% COVDATA.COVIDs = [COVDATA.COVIDs;addCOVDATA.COVIDs(index)];
% COVDATA.COVNames = [COVDATA.COVNames;addCOVDATA.COVNames(index)];
% COVDATA.DATA = [COVDATA.DATA addCOVDATA.DATA(ind12,index)];
% % adding second set of covariates
% addCOVDATA = COVDATA2;
% [ind12,~] = vlookupFast(COVDATA.IID,addCOVDATA.IMGID);
% index = find(~ismember(addCOVDATA.COVIDs,COVDATA.COVIDs));
% COVDATA.COVIDs = [COVDATA.COVIDs;addCOVDATA.COVIDs(index)];
% COVDATA.COVNames = [COVDATA.COVNames;addCOVDATA.COVNames(index)];
% COVDATA.DATA = [COVDATA.DATA addCOVDATA.DATA(ind12,index)];
% COVDATA.nCOV = length(COVDATA.COVIDs);
% save('COVDATA123','COVDATA','-v7.3');
% %% CONVERSIONS & SELECTIONS
% load('COVDATA123');
% % Subjects info
% COV.IID = COVDATA.IID;
% COV.FID = COVDATA.FID;
% % Genetic PCs
% Index = 1:20;%index of genetic PCs
% COV.Names = COVDATA.COVNames(Index);
% COV.IDs = COVDATA.COVIDs(Index);
% Data = cell2mat(COVDATA.DATA(:,Index));
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = Data;
% % Genetic Sex
% Index = 27;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% Data = COVDATA.DATA(:,Index);
% Data = listofstr2num(Data);
% Data = cell2mat(Data);
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% % Age
% Index = 26;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% Data = COVDATA.DATA(:,Index);
% Data = listofstr2num(Data);
% Data = cell2mat(Data);
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% % Age Squared
% COV.Names = [COV.Names; {[COVDATA.COVNames{Index} ' Squared']}];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% COV.DATA = [COV.DATA Data.^2];
% % Height
% Index = 32;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% Data = COVDATA.DATA(:,Index);
% Data = listofstr2num(Data);
% Data = cell2mat(Data);
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% % Weight
% Index = 30;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% Data = COVDATA.DATA(:,Index);
% Data = listofstr2num(Data);
% Data = cell2mat(Data);
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% % DIASTOLIC BLOOD PRESSURE
% Index = 35;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% Data = COVDATA.DATA(:,Index);
% Data = listofstr2num(Data);
% Data = cell2mat(Data);
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% % SYSTOLIC BLOOD PRESSURE
% Index = 36;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% Data = COVDATA.DATA(:,Index);
% Data = listofstr2num(Data);
% Data = cell2mat(Data);
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% %
% COV.GENOIndex = 1:length(COV.Names);
% %
% % VOLUMETRIC SCALING
% Index = 40;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% Data = COVDATA.DATA(:,Index);
% Data = listofstr2num(Data);
% Data = cell2mat(Data);
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% % POSITION
% IndexList = 41:44;
% for i=1:1:length(IndexList)
%     Index = IndexList(i);
%     COV.Names = [COV.Names; COVDATA.COVNames(Index)];
%     COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
%     Data = COVDATA.DATA(:,Index);
%     Data = listofstr2num(Data);
%     Data = cell2mat(Data);
%     disp(COVDATA.COVNames{Index});
%     disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
%     Data = CenterDataInMatrix(Data,true);
%     COV.DATA = [COV.DATA Data];
% end
% % ASSESMENT TIME
% Index = 24;
% COV.Names = [COV.Names; COVDATA.COVNames(Index)];
% COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
% TMP = COVDATA.DATA(:,Index);
% Data = nan*zeros(length(TMP),1);
% for i=1:1:length(TMP)
%     str = TMP{i};
%     if isempty(str), continue; end
%     ind = strfind(str,'-');
%     Y = str2num(str(1:ind(1)-1)); %#ok<*ST2NM>
%     M = str2num(str(ind(1)+1:ind(2)-1));
%     D = str2num(str(ind(2)+1:end));
%     Data(i) = datenum(Y,M,D);
% end
% disp(COVDATA.COVNames{Index});
% disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
% Data = CenterDataInMatrix(Data,true);
% COV.DATA = [COV.DATA Data];
% % ASSESMENT CENTER
% Index = 23;
% name = COVDATA.COVNames(Index);
% TMP = COVDATA.DATA(:,Index);
% TMP = listofstr2num(TMP);
% TMP = cell2mat(TMP);
% ass = unique(TMP);
% for i=1:1:length(ass)-1% all but last (to avoid colinearity)
%     %i=1;
%     COV.Names = [COV.Names; {[COVDATA.COVNames{Index} ' ' num2str(ass(i))]}];
%     COV.IDs = [COV.IDs; {[COVDATA.COVIDs{Index} ' ' num2str(ass(i))]}];
%     Data = double(TMP==ass(i));
%     COV.DATA = [COV.DATA Data];
% end
% %
% COV.PHENOIndex = 1:length(COV.Names);
% %
% save('COVDATA2USE','COV','-v7.3');
% %% CHECK FOR OUTLIERS
% load('COVDATA2USE');
% [nSUB,nVAR] = size(COV.DATA);
% outliers = zeros(nSUB,nVAR);
% nrOutliers = zeros(1,nVAR);
% for i=1:1:nVAR
%     %i=24;
%     val = COV.DATA(:,i);
%     TF = isoutlier(val,'mean','ThresholdFactor',6);
%     index = find(TF==1);
%     nrOutliers(i) = length(index);
%     outliers(:,i) = TF(:);
% end
% COV.outliers = outliers;
% COV.nrOutliers = nrOutliers;
% %save('COVDATA2USE','COV','-v7.3');
% % There are like 11 people with an excessive wheith on 20.000, but I think
% % this is fine, and not unlikely, there values are large but not impossible
% %% REMOVING OUTLIERS
% load('COVDATA2USE');
% % first print
% for i=1:length(COV.IDs)
%     disp(num2str(i));
%     disp(COV.Names{i});
%     disp(num2str(COV.nrOutliers(i)));
% end
% % I would only eliminate te weight people, but then again this has not the
% % same effect on brain structure as it has on facial structure
% OUT = sum(COV.outliers(:,25:33),2);
% INLIERS = find(OUT==0);
% COV.IndexInliers = INLIERS;
% %save('COVDATA2USE','COV','-v7.3');
% COV.IID = COV.IID(INLIERS);
% COV.FID = COV.FID(INLIERS);
% COV.DATA = COV.DATA(INLIERS,:);
% COV = rmfield(COV,'outliers');
% COV = rmfield(COV,'nrOutliers');
% COV = rmfield(COV,'IndexInliers');
% %save('COVDATAINLIERS','COV','-v7.3');
%
