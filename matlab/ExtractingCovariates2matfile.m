%%
close all;clear all;
restoredefaultpath;
addpath(genpath('AIDFUNCTIONS'));

DATA_DIR = '../SAMPLE_DATA/';
PHENO_DIR = [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/PHENOTYPES/'];
GENO_DIR = [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/GENOTYPES/'];

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
%%
csvpath0 = fullfile(primarycovpath, 'ukb25183.csv');
[parsedData0, ids0, cols0] = readcsv(csvpath0, ids, COVIDs);
%%

csvpath1 = fullfile(primarycovpath, 'ukb35648.csv');
[parsedData1, ids1, cols1] = readcsv(csvpath1, ids, COVIDs);
%%
csvpath2 = fullfile(secondarycovpath, 'ukb44263', 'ukb44263.csv');
[parsedData2, ids2, cols2] = readcsv(csvpath2, ids, COVIDs);
%% TODO, adapt the following to connect the two datasets, check ids first!
[finalData, finalIds, finalCols] = mergeFiles(parsedData0, ids0, cols0, parsedData1, ids1, cols1);
[finalData, finalIds, finalCols] = mergeFiles(finalData, finalIds, finalCols, parsedData2, ids2, cols2);
%% Apply columns reduction again, if sister columns were added in the stages above (belonging to the same id)
[filteredCols, fCOVIDsInds, colInds]  = findCovIds(finalCols, COVIDs);
finalData = reduce(finalData, colInds);

COVDATA.nCOV = length(filteredCols);
COVDATA.COVIDs = filteredCols;
COVDATA.COVNames = CovNames(cell2mat(fCOVIDsInds));

COVDATA.nSUB = length(finalIds);
COVDATA.IMGID = finalIds;
COVDATA.DATA = finalData;


%% IMPORTING GENETIC AXIS
in = load(fullfile(GENO_DIR, 'UKB_EUR_16875_PCs_2.mat'));
nPC = size(in.AA,2);
pcIID = in.IID;
pcCOVIDs = cell(nPC,1);
pcCOVNames = cell(nPC,1);
for i=1:nPC
    pcCOVIDs{i} = ['00' num2str(i)];
    pcCOVNames{i} = ['KUL GENO PC ' num2str(i)];
end
find(pc)


COVDATA3.DATA = num2cell(in.AA);

[ind12,~] = vlookupFast(COVDATA.IID,addCOVDATA.IMGID);
index = find(~ismember(addCOVDATA.COVIDs,COVDATA.COVIDs));
COVDATA.COVIDs = [COVDATA.COVIDs;addCOVDATA.COVIDs(index)];
COVDATA.COVNames = [COVDATA.COVNames;addCOVDATA.COVNames(index)];
COVDATA.DATA = [COVDATA.DATA addCOVDATA.DATA(ind12,index)];


function [finalData, IID, COLS] = mergeFiles(oldData, oldIds, oldCols, newData,newIds, newCols)
[ids1Ind, ids2Ind ] = find(str2double(oldIds) == str2double(newIds)');

extraOldRows = 1:length(oldIds);
extraOldRows = extraOldRows(~ismember(extraOldRows,ids1Ind));
extraNewRows = 1:length(newIds);
extraNewRows = extraNewRows(~ismember(extraNewRows,ids2Ind));

idsInd = [ids1Ind; extraOldRows'; extraNewRows'];
IID = [oldIds(ids1Ind); oldIds(extraOldRows); newIds(extraNewRows)];
[cols1Ind, cols2Ind] = find(oldCols == newCols');
extraOldCols = 1:length(oldCols);
extraOldCols = extraOldCols(~ismember(extraOldCols,cols1Ind));
extraNewCols = 1:length(newCols);
extraNewCols = extraNewCols(~ismember(extraNewCols,cols2Ind));
% colsInds = [cols1Ind; extraOldCols'; extraNewCols'];
COLS = [oldCols(cols1Ind); oldCols(extraOldCols); newCols(extraNewCols)];
oldData2 = oldData(ids1Ind, cols1Ind);
newData2 = newData(ids2Ind, cols2Ind);

upd_mask = ~cellfun(@(c) (isscalar(c) && isnan(c)) || isempty(c), newData2);
finalData = oldData2;
finalData(upd_mask) = newData2(upd_mask);

finalData = [finalData, oldData(ids1Ind, extraOldCols), newData(ids2Ind, extraNewCols)];
if ~isempty(extraOldRows)
    remOld =  zeros(length(extraOldRows), length(COLS));
    remOld(:) = nan;
    remOld = num2cell(remOld);
    remOld(:, 1:(length(COLS)-length(extraNewCols))) = oldData(extraOldRows, [cols1Ind; extraOldCols']);
    finalData = [finalData; remOld];
end
if ~isempty(extraNewRows)
    remNew = zeros(length(extraNewRows), length(COLS));
    remNew(:) = nan;
    remNew = num2cell(remNew);
    remNew(:,  [1:(length(cols1Ind)), (length(COLS)-length(extraNewCols) + 1) : length(COLS) ]) = ...
        newData(extraNewRows, [cols2Ind; extraNewCols']);
    finalData = [finalData; remNew];
end
end



function  [filteredCols, fCOVIDsInds, colInds] = findCovIds(headers, COVIDs)
    if any(arrayfun(@isstring, headers))
        colIndices = arrayfun(@(x) find(startsWith(headers, ['"', num2str(x) '-'])), COVIDs , 'UniformOutput', false);
    else
        colIndices = arrayfun(@(x) find(headers==x), COVIDs , 'UniformOutput', false);
    end
    found = ~cellfun(@isempty, colIndices);
    fCOVIDs = COVIDs(found);
    
    filteredCols = fCOVIDs;
    fCOVIDsInds = colIndices(found);
    num = 1:length(fCOVIDsInds);
    colInds = cellfun(@(x, ind)(repmat(ind, length(x), 1)), fCOVIDsInds, num2cell(num)', 'UniformOutput', false);
    colInds = vertcat(colInds{:}); % No need to add first column index, as it is removed after from data
end



function [parsedData, filteredIds, filteredCols ]= readcsv(csvpath, ids, COVIDs, debugNLines)
    fid = fopen(csvpath);
    headers = strsplit(fgetl(fid), ',');
    numCol = length( strsplit(fgetl(fid), ','));
    fclose(fid);
    [filteredCols, fCOVIDsInds, colInds] = findCovIds(headers, COVIDs);
    fmt_parts = repmat({'%*s'}, 1, numCol);
    fmt_parts(horzcat(fCOVIDsInds{:})) = {'%q'};
    fmt_parts(1) = {'%q'};
    fmt = strjoin(fmt_parts, '');
    fid = fopen(csvpath);
    if nargin == 3
        data = textscan(fid, fmt, 'Delimiter', ',', 'Headerlines', 1);
    else
        data = textscan(fid, fmt, debugNLines, 'Delimiter', ',', 'Headerlines', 1);
    end
    data{1} = num2cell(str2double(data{1}));
    fclose(fid);

    data = [data{:}];

    all_ids =[data{:, 1}];
    [all_ids, sortInds] = sort(all_ids);
    data = data(sortInds, 2:end);
    [to_keep, mapped_to] = find(str2double(ids) == all_ids);
    filteredIds = ids(to_keep);
    data = data(mapped_to, :);
    parsedData = reduce(data, colInds);
end


%%
function parsedData = reduce(data, colInds)
    uniqueInds = unique(colInds);
    parsedData = zeros(size(data, 1), length(uniqueInds));
    parsedData(:) = nan;
    parsedData = num2cell(parsedData);
    for c = 1:length(uniqueInds)
        subset = data(:, colInds == c);
        ret = zeros(size(data, 1), 1);
        ret(:) = nan;
        ret = num2cell(ret);
        for col = 1:size(subset,2)
            mask = ~cellfun(@isempty, subset(:,col));
            ret(mask) = subset(mask, col);
        end
        parsedData(:, c) = ret;
    end
end


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
