function outFiles=runCCA(pheno, jobsFiles, outFiles, incMem)
%pheno: path to the phenotype
% jobsFiles: a list of input jobs files, or their directory
% outFiles: if provided, a list of output files or the output directory
% incMem : increase memory usage for the sake of speed, defaults to false,
% as it does not work well when the dataset is not imputed
if nargin<4
    incMem=0;
end
if isa(pheno,"char")
    pheno=load(pheno).PHENO;
end
if isa(jobsFiles,'char')
    jobsFiles = dir(jobsFiles);
    inpfolder = jobsFiles(1).folder;
    numFiles =sum(cellfun(@(x)(endsWith(x,'.mat')), {jobsFiles.name}));
    names = strcat(strsplit(num2str(1:numFiles)) ,'.mat');
    jobsFiles =   strcat([inpfolder '/'], names);
else
    parts = strsplit(jobsFiles{1}, filesep);
    inpfolder = strjoin(parts(1:end-1), filesep);
    names = cell(length(jobsFiles));
    for i=1:length(jobsFiles)
        path = strsplit(jobsFiles{i}, filesep);
        names{i} = path(end);
    end
end

if nargin < 3
    outfolder = [inpfolder '_out/'];
else
    if isa(outFiles, 'char')
        outfolder = outFiles;
    end
end

if nargin <3 || isa(outFiles,'char')
    if ~isfolder(outfolder), mkdir(outfolder); end
    outFiles = strcat(outfolder, names);
end

phenoPartNum = length(pheno);
Q2Full = nan;
T22Full = nan;
rankYFull = nan;
if incMem
    Q2Full = cell(length(pheno),1);
    T22Full = cell(length(pheno),1);
    rankYFull = cell(length(pheno),1);
    parfor i=1:length(pheno)
        [~, Q2Full{i}, T22Full{i}, rankYFull{i}] = vl_mycanoncorr([], pheno{i});
    end
end
ME = [];
fallback = ~exist('pagesvd','builtin') || ~exist('pagemtimes', 'builtin');
try
    ppb = ParforProgressbar(length(jobsFiles),  'showWorkerProgress', true);
    parfor blockCnt=1: length(jobsFiles)
        inp = load(jobsFiles{blockCnt});
        genoBlock = inp.genoBlock;
        blockIntervals = inp.blockIntervals;
        blockSizes = blockIntervals(:,2) - blockIntervals(:,1) + 1;
        gId =  zeros(sum(blockSizes),1);
        gId(blockIntervals(:, 1) - blockIntervals(1,2) + 1) = 1;
        gId = cumsum(gId);
        genoBlock = splitapply( @(x){x'}, genoBlock', gId );
        currBlockSize = size(genoBlock,1);
        finite_mask = cellfun(@(x)all(x~=255,2), genoBlock,'UniformOutput',false);
        finite_mask = cat(2,finite_mask{:})';
        [patterns, ~, patternsInds] = unique(finite_mask,'rows');
        intChiSqSignificancePattern = ones(currBlockSize, phenoPartNum);
        intChiSqPattern = zeros(currBlockSize,phenoPartNum);
        intDfPattern = zeros(currBlockSize,phenoPartNum);
        for pcnt=1:size(patterns,1)
            pattern = patterns(pcnt,:);
            patternFlag= patternsInds == pcnt;
            patternSizes = blockSizes(patternFlag);
            patternNSnps = sum(patternFlag);
            nomissing = all(patternFlag);
            if ~nomissing
                genoPattern =cellfun(@(x)x(pattern,:), genoBlock(patternFlag), 'UniformOutput', false);
            else
                genoPattern = genoBlock;
            end
            for phenoPartInd=1:length(pheno)
                phenoPart = pheno{phenoPartInd};
                phenoPart = phenoPart(pattern, :);
                if patternNSnps==1
                    if size(genoPattern{1},1) ~= size(phenoPart,1)
                        disp('Error')
                    end
                    st = vl_mycanoncorr(double(genoPattern{1}), phenoPart);
                    intChiSqPattern(patternFlag, phenoPartInd)=  st.chisq(1);
                    intChiSqSignificancePattern(patternFlag, phenoPartInd) = st.pChisq(1);
                    intDfPattern(patternFlag, phenoPartInd) = st.df1(1);
                    continue
                end
                if ~nomissing || ~incMem
                    [~, Q2, T22, rankY] = vl_mycanoncorr([], phenoPart);
                else
                    Q2 = Q2Full{phenoPartInd};
                    T22 = T22Full{phenoPartInd};
                    rankY = rankYFull{phenoPartInd};
                end
                if fallback
                    [intChiSqPattern(patternFlag, phenoPartInd), ...
                        intChiSqSignificancePattern(patternFlag, phenoPartInd) , ...
                        intDfPattern(patternFlag, phenoPartInd)] = fallbackFunc(genoPattern, phenoPart, Q2, T22, rankY);
                    continue
                end
                intChiSqSignificanceSize = ones(patternNSnps, 1);
                intChiSqSize =  zeros(patternNSnps, 1);
                dfsSize =  zeros(patternNSnps,1);
                for s=1:max(patternSizes)
                    sizeFlag = (patternSizes== s);
                    if ~any(sizeFlag)
                        continue
                    end
                    part = genoPattern(sizeFlag);
                    Xs = double(cat(3,part{:}));
                    try
                        [intChiSqSize(sizeFlag),  intChiSqSignificanceSize(sizeFlag), dfsSize(sizeFlag)] = vl_ccachisq1(Xs, nan,Q2, T22, rankY);
                    catch Exception
                        disp(Exception);
                    end
                end
                intChiSqSignificancePattern(patternFlag, phenoPartInd) = intChiSqSignificanceSize;
                intChiSqPattern(patternFlag, phenoPartInd) = intChiSqSize;
                intDfPattern(patternFlag, phenoPartInd) = dfsSize;
            end
        end
        saveBlockFunc(outFiles{blockCnt}, intChiSqSignificancePattern, intChiSqPattern, intDfPattern);
        ppb.increment();
    end
catch ME
end
delete(ppb);
if ~isempty(ME)
    rethrow(ME);
end
end


function saveBlockFunc(outfile, chisqSignificance, chisq,df)
save(outfile, 'chisqSignificance','chisq','df')
end


function [intChiSq, intChiSqSignificance, dfs] = fallbackFunc(geno, phenoPart, Q2, T22, rankY)
%If not supporting parallel SVD
ngeno = length(geno);
intChiSq = zeros(ngeno,1);
intChiSqSignificance =  ones(ngeno,1);
dfs = zeros(ngeno, 1);
for i=1:length(geno)
    genoPart = double(geno{i});
    st = vl_mycanoncorr(genoPart, phenoPart, Q2, T22, rankY);
    intChiSq(i) =  st.chisq(1);
    intChiSqSignificance(i) = st.pChisq(1);
    dfs(i) = st.df1(1);
end
end

