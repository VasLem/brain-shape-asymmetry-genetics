%%
close all;clear all;
restoredefaultpath;
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('/opt/SNPLIB/'));
rmpath('/opt/SNPLIB/mexfiles/');% to remove the functions that are causing matlab to crash
addpath(genpath('SNPLIB-master/mexfiles/'))% where I stored the re-mexed files
DATA_DIR = '../SAMPLE_DATA/';
PHENO_DIR = [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/PHENOTYPES/'];
GENO_DIR = [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/GENOTYPES/PLINK/'];
COV_DIR  =  [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/COVARIATES'];
if ~isfolder(COV_DIR), mkdir(COV_DIR); end
THREADS = 8;
primarycovpath = '/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/COVARIATES';
secondarycovpath = '../SAMPLE_DATA/IMAGEN/BRAIN/BATCH2_DATA/';
% covpath = [datapath 'COVARIATES/'];
% phenopath = [datapath 'PHENOTYPING/'];
% load([phenopath 'FinalImageIDs_20407subjs']);%loads imageIDs
% IMGID = num2cell(imageIDs);
% IMGID = listofnum2str(IMGID);
% cd(covpath);

%% GETTING RELEVANT DATA ENTRY FIELDS AND IIDS
x = readtable(fullfile(primarycovpath,'Final_list_covariates_GWAS.xlsx' ));
COVIDs = x.Var1;
CovNames = x.Var2;
nCOV = length(COVIDs);

ids = load(fullfile(PHENO_DIR, 'BATCH2_2021_DATA_IID')).iids;
%% READ AND REDUCE CSV FILES
% csvpath0 = fullfile(primarycovpath, 'ukb25183.csv');
% [parsedData0, ids0, cols0] = readcsv(csvpath0, ids, COVIDs);


% csvpath1 = fullfile(primarycovpath, 'ukb35648.csv');
% [parsedData1, ids1, cols1] = readcsv(csvpath1, ids, COVIDs);
%%
disp("Reading 1st csv file..")
csvpath1 = fullfile(secondarycovpath, 'ukb44267', 'ukb44267.csv');
[parsedData1, ids1, cols1] = readcsv(csvpath1, ids, COVIDs);
%%
disp("Reading 2nd csv file..")
csvpath2 = fullfile(secondarycovpath, 'ukb44263', 'ukb44263.csv');
[parsedData2, ids2, cols2] = readcsv(csvpath2, ids, COVIDs);
%% MERGE REDUCED FILES
disp("Merging files..")
[finalData, finalIds, finalCols] = mergeFiles(parsedData1, ids1, cols1, parsedData2, ids2, cols2);
% [finalData, finalIds, finalCols] = mergeFiles(parsedData0, ids0, cols0, finalData, finalIds, finalCols);
%% SORT COLUMNS AND APPLY COLUMN REDUCTION AGAIN
% if sister columns were added in the stages above (belonging to the same id)
% For example if 50-1 resided in first csv and 50-2 in the second, merge
% them into 50
[filteredCols, fCOVIDsInds, ~]  = findCovIds(finalCols, COVIDs);
finalData = reduce(finalData, fCOVIDsInds);
%% CREATE FIRST PART OF THE FINAL STRUCTURE

COVDATA.nCOV = length(filteredCols);
COVDATA.COVIDs = cellstr(num2str(filteredCols));
COVDATA.COVNames = convertCharsToStrings(x(ismember(x.Var1, filteredCols), :).Var2);

COVDATA.nSUB = length(finalIds);
COVDATA.IMGID = vertcat(finalIds{:});
COVDATA.DATA = finalData;
T = array2table(finalData, "VariableNames", cellfun(@(x) x(1:min(30, length(x)) ),  convertCharsToStrings(x(ismember(x.Var1, filteredCols), :).Var2), 'UniformOutput', 0));
summary(T);

%% IMPORTING PCs


 % Use small chr 21, to get iids from.

obj = SNPLIB();
obj.nThreads = THREADS;
[snps, samples] = obj.importPLINKDATA([GENO_DIR, 'ukb_img_maf0.01_geno0.5_hwe1e-6_sel16875_rmrel_chr21']);
pcIID = samples.IID;

% Load eigen vectors (PCs)

AA = load(fullfile(GENO_DIR, 'ukb_img_maf0.01_geno0.5_hwe1e-6_sel16875_rmrel_pca.eigenvec'));
%%
[genoIndIds, covIndIds] = find(str2double(pcIID) == str2double(COVDATA.IMGID)');
nPC = size(AA,2);

pcCOVIDs = cell(nPC,1);
pcCOVNames = cell(nPC,1);
for i=1:nPC
    pcCOVIDs{i} = ['00' num2str(i)];
    pcCOVNames{i} = ['KUL GENO PC ' num2str(i)];
end
COVDATA.COVNames = [pcCOVNames; COVDATA.COVNames ];
COVDATA.COVIDs = [pcCOVIDs; COVDATA.COVIDs];
COVDATA.DATA = [AA(genoIndIds, :) COVDATA.DATA(covIndIds, :)];
COVDATA.IMGID = pcIID(genoIndIds) ;

%% Genetic PCs
Index = 1:20;%index of genetic PCs
COV.Names = COVDATA.COVNames(Index);
COV.IDs = COVDATA.COVIDs(Index);
COV.IID = COVDATA.IMGID;
Data = str2double(COVDATA.DATA(:,Index));
Data = CenterDataInMatrix(Data,true);
COV.DATA = Data;
%% Genetic Sex

%%
Index = find(COVDATA.COVNames  == 'Genetic sex');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
Data = COVDATA.DATA(:,Index);
Data = str2double(Data);
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%% Age
Index =  find(COVDATA.COVNames  == 'Age when attended assessment centre');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
Data = COVDATA.DATA(:,Index);
Data = str2double(Data);
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%% Age Squared
COV.Names = [COV.Names; {[COVDATA.COVNames{Index} ' Squared']}];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
COV.DATA = [COV.DATA Data.^2];
%% Height
Index =  find(COVDATA.COVNames  == 'Standing height');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
Data = COVDATA.DATA(:,Index);
Data = str2double(Data);
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%% Weight
Index =  find(COVDATA.COVNames  == 'Weight');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
Data = COVDATA.DATA(:,Index);
Data = str2double(Data);
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%% DIASTOLIC BLOOD PRESSURE
Index =  find(COVDATA.COVNames  == 'Diastolic blood pressure, automated reading');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
Data = COVDATA.DATA(:,Index);
Data = str2double(Data);
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%% SYSTOLIC BLOOD PRESSURE
Index =  find(COVDATA.COVNames  == 'Systolic blood pressure, automated reading');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
Data = COVDATA.DATA(:,Index);
Data = str2double(Data);
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%
COV.GENOIndex = 1:length(COV.Names);
%
%% VOLUMETRIC SCALING
Index =  find(COVDATA.COVNames  == 'Volumetric scaling from T1 head image to standard space');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
Data = COVDATA.DATA(:,Index);
Data = str2double(Data);
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%% POSITION
IndexList =  find(ismember( COVDATA.COVNames, ["X-position of centre-of-gravity of brain mask in scanner co-ordinates",...
    "Y-position of back of brain mask in scanner co-ordinates",...
    "Z-position of centre-of-gravity of brain mask in scanner co-ordinates",...
    "Z-position of table/coil in scanner co-ordinates"]));

for i=1:1:length(IndexList)
    Index = IndexList(i);
    COV.Names = [COV.Names; COVDATA.COVNames(Index)];
    COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
    Data = COVDATA.DATA(:,Index);
    Data = str2double(Data);
    disp(COVDATA.COVNames{Index});
    disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
    Data = CenterDataInMatrix(Data,true);
    COV.DATA = [COV.DATA Data];
end
%% ASSESMENT TIME
Index = find(COVDATA.COVNames  == 'Date of attending Assessment centre');
COV.Names = [COV.Names; COVDATA.COVNames(Index)];
COV.IDs = [COV.IDs; COVDATA.COVIDs(Index)];
TMP = COVDATA.DATA(:,Index);
Data = nan*zeros(length(TMP),1);
for i=1:1:length(TMP)
    str = TMP{i};
    if isempty(str), continue; end
    ind = strfind(str,'-');
    Y = str2num(str(1:ind(1)-1)); %#ok<*ST2NM>
    M = str2num(str(ind(1)+1:ind(2)-1));
    D = str2num(str(ind(2)+1:end));
    Data(i) = datenum(Y,M,D);
end
disp(COVDATA.COVNames{Index});
disp(['Missing Data: ' num2str(length(find(isnan(Data))))]);
Data = CenterDataInMatrix(Data,true);
COV.DATA = [COV.DATA Data];
%% ASSESMENT CENTER
Index = find(COVDATA.COVNames  =='UK Biobank assessment centre');
name = COVDATA.COVNames(Index);
TMP = COVDATA.DATA(:,Index);
TMP = str2double(TMP);
ass = unique(TMP);
for i=1:1:length(ass)-1% all but last (to avoid colinearity)
    %i=1;
    COV.Names = [COV.Names; {[COVDATA.COVNames{Index} ' ' num2str(ass(i))]}];
    COV.IDs = [COV.IDs; {[COVDATA.COVIDs{Index} ' ' num2str(ass(i))]}];
    Data = double(TMP==ass(i));
    COV.DATA = [COV.DATA Data];
end
%
COV.PHENOIndex = 1:length(COV.Names);
%
save(fullfile(COV_DIR,  'COVDATA2USE'),'COV','-v7.3');
%% CHECK FOR OUTLIERS
load('COVDATA2USE');
[nSUB,nVAR] = size(COV.DATA);
outliers = zeros(nSUB,nVAR);
nrOutliers = zeros(1,nVAR);
for i=1:1:nVAR
    %i=24;
    val = COV.DATA(:,i);
    TF = isoutlier(val,'mean','ThresholdFactor',6);
    index = find(TF==1);
    nrOutliers(i) = length(index);
    outliers(:,i) = TF(:);
end
COV.outliers = outliers;
COV.nrOutliers = nrOutliers;
%save('COVDATA2USE','COV','-v7.3');
% There are like 11 people with an excessive wheith on 20.000, but I think
% this is fine, and not unlikely, there values are large but not impossible
%% REMOVING OUTLIERS
% first print
for i=1:length(COV.IDs)
    disp(num2str(i));
    disp(COV.Names{i});
    disp(num2str(COV.nrOutliers(i)));
end
% I would only eliminate te weight people, but then again this has not the
% same effect on brain structure as it has on facial structure
%%
OUT = sum(COV.outliers(:,25:33),2);
INLIERS = find(OUT==0);
COV.IndexInliers = INLIERS;
%save('COVDATA2USE','COV','-v7.3');
COV.IID = COV.IID(INLIERS);
COV.DATA = COV.DATA(INLIERS,:);
COV = rmfield(COV,'outliers');
COV = rmfield(COV,'nrOutliers');
COV = rmfield(COV,'IndexInliers');
save(fullfile(COV_DIR, 'COVDATAINLIERS'),'COV','-v7.3');
%%

function [finalData, IID, COLS] = mergeFiles(oldData, oldIds, oldCols, newData,newIds, newCols)
[ids1Ind, ids2Ind ] = find(str2double([oldIds{:}]) == str2double([newIds{:}])');

extraOldRows = 1:length(oldIds);
extraOldRows = extraOldRows(~ismember(extraOldRows,ids1Ind));
extraNewRows = 1:length(newIds);
extraNewRows = extraNewRows(~ismember(extraNewRows,ids2Ind));

% idsInd = [ids1Ind; extraOldRows'; extraNewRows'];
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

upd_mask = newData2~='';
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
    remNew = strings(length(extraNewRows), length(COLS));
    remNew(:,  [1:(length(cols1Ind)), (length(COLS)-length(extraNewCols) + 1) : length(COLS) ]) = ...
        newData(extraNewRows, [cols2Ind; extraNewCols']);
    finalData = [finalData; remNew];
end
end



function  [filteredCols, fCOVIDsInds, fCOVIDDenseInds] = findCovIds(headers, COVIDs)
    if  iscell(headers)
        colIndices = arrayfun(@(x) find(startsWith(headers, ['"', num2str(x) '-'])), COVIDs , 'UniformOutput', false);
    else
        colIndices = arrayfun(@(x) find(headers==x), COVIDs , 'UniformOutput', false);
    end
    found = ~cellfun(@isempty, colIndices);
    fCOVIDs = COVIDs(found);
    
    filteredCols = fCOVIDs;
    fCOVIDsInds = colIndices(found);
    colInds = vertcat([fCOVIDsInds{:}]);
    [~, a] = sort(colInds);
    [~, denseInds] = sort(a);
    fCOVIDDenseInds = cell(length(fCOVIDsInds), 1);
    c = 1;
    for c1 = 1: length(fCOVIDsInds)
        for c2 = 1:length(fCOVIDsInds{c1})
            fCOVIDDenseInds{c1}(c2) = denseInds(c);
            c = c + 1;
        end
    end

end



function [parsedData, filteredIds, filteredCols ]= readcsv(csvpath, ids, COVIDs, debugNLines)
    outname = tempname;
    disp("Adding tabs delimiters to original file..")
    system(sprintf("sed 's/%s/%s/g' %s > %s", '","', '"\t"', csvpath, outname));
    fid = fopen(outname);
    headers = strsplit(fgetl(fid), '\t');
    numCol = length(headers);
    fprintf("Number of headers in file:%d\n", numCol);
    fclose(fid);
    disp("Finding covariates Ids groups inside file (eg 50-1, 50-2 columns refer to ID 50) ..")
    [filteredCols, fCOVIDsInds, denseInds] = findCovIds(headers, COVIDs);
    fprintf("%d out of %d covariates were found\n", length(filteredCols), length(COVIDs))
    disp("Constructing formatter string to scan file..")
    fmt_parts = repmat({'%*s'}, 1, numCol);
    fmt_parts(horzcat(fCOVIDsInds{:})) = {'%q'};
    fmt_parts(1) = {'%q'};
    fmt = strjoin(fmt_parts, '');
    fid = fopen(outname);
    disp("Getting found columns from file..")
    if nargin == 3
        data = textscan(fid, fmt, 'Delimiter', '\t', 'Headerlines', 1);
    else
        data = textscan(fid, fmt, debugNLines,  'Delimiter',  '\t' , 'Headerlines', 1);
    end
    fclose(fid);
    delete(outname);
    all_ids = data{1};
    [~, mapped_to] = find(str2double(ids) == str2double(all_ids)');
    data = [data{:}];

    
    filteredIds = all_ids(mapped_to);
    fprintf("Found %d out of %d samples ids\n", length(mapped_to), length(ids));
    [filteredIds, sortInds] = sort(filteredIds);
    filteredIds = num2cell(filteredIds);
    data = data(sortInds, 2:end);
    data = convertCharsToStrings(data);
    disp("Aggregating ID groups, to retain newest non missing information..")
    parsedData = reduce(data, denseInds);
end

function parsedData = reduce(data, fCOVIDsInds)
    parsedData = strings(size(data, 1), length(fCOVIDsInds));
    for c = 1:length(fCOVIDsInds)
        subset = data(:, fCOVIDsInds{c});
        for col = 1:size(subset,2)
            mask = subset(:,col) ~= "";
            parsedData(mask, c) = subset(mask, col);
        end
    end
end

