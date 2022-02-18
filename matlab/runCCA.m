function [stats, intStats] = runCCA(pheno, geno, intervals, intIdVec, blockSize)
intChiSqSignificance = zeros(size(intervals, 1),1);
intChiSq = zeros(size(intervals, 1),1);
intDf = zeros(size(intervals,1),1);
phenoPartNum = length(pheno);
if nargin < 4
    disp("Computing intIdVec (expanded version of intervals)..")
    intIdVec = zeros(intervals(end,2),1);
    intIdVec(intervals(:,1)) = 1;
    intIdVec = cumsum(intIdVec);
end
if ~iscell(geno)
    disp("Converting geno to cell array..")
    geno = splitapply( @(x){x'}, geno', intIdVec );
end

genoSizes = intervals(:,2) - intervals(:,1) + 1;

total = size(geno,1);
blocksN = ceil(total/blockSize);
intChiSqSignificanceBlock = cell(blocksN,1);
intChiSqBlock =  cell(blocksN,1);
intDfBlock =  cell(blocksN,1);
genoBlocks = cell(blocksN, 1);
minInds = zeros(blockSize,1);
maxInds = zeros(blockSize,1);
for blockCnt=1: blocksN
    minInds(blockCnt) = (1 + (blockCnt - 1) * blockSize);
    maxInds(blockCnt) = min(total, (blockCnt * blockSize));
    genoBlocks{blockCnt} = geno(minInds(blockCnt):maxInds(blockCnt));
end
Q2Full = cell(length(pheno),1);
T22Full = cell(length(pheno),1);
rankYFull = cell(length(pheno),1);
parfor i=1:length(pheno)
    [~, Q2Full{i}, T22Full{i}, rankYFull{i}] = vl_mycanoncorr([], pheno{i});
end
ME = [];
try
    %%
    h = waitbar(0,'Computing CCA...');
    clk = clock;
    for blockCnt=1: blocksN
        currBlockSize = maxInds(blockCnt) - minInds(blockCnt) + 1;
        genoBlock = genoBlocks{blockCnt};
        finite_mask = cellfun(@(x)all(x~=255,2), genoBlock,'UniformOutput',false);
        finite_mask = cat(2,finite_mask{:})';
        [patterns, ~, patternsInds] = unique(finite_mask,'rows');
        fallback = ~exist('pagesvd','builtin') || ~exist('pagemtimes', 'builtin');
        blockSizes = genoSizes(minInds(blockCnt):maxInds(blockCnt));
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
            parfor phenoPartInd=1:length(pheno)
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
                if ~nomissing
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
                    [intChiSqSize(sizeFlag),  intChiSqSignificanceSize(sizeFlag), dfsSize(sizeFlag)] = vl_ccachisq1(Xs, nan,Q2, T22, rankY);
                end
                intChiSqSignificancePattern(patternFlag, phenoPartInd) = intChiSqSignificanceSize;
                intChiSqPattern(patternFlag, phenoPartInd) = intChiSqSize;
                intDfPattern(patternFlag, phenoPartInd) = dfsSize;
            end
        end
        intChiSqSignificanceBlock{blockCnt} = intChiSqSignificancePattern;
        intChiSqBlock{blockCnt} = intChiSqPattern;
        intDfBlock{blockCnt} = intDfPattern;
        if blockCnt ==1
            is = etime(clock,clk);
            esttime = is * blocksN;
        end
        h = waitbar(blockCnt/blocksN,h,...
            ['remaining time =',num2str(esttime-etime(clock,clk),'%4.1f'),'sec' ]);
    end
    intChiSqSignificance = cat(1, intChiSqSignificanceBlock{:});
    intChiSq = cat(1, intChiSqBlock{:});
    intDf = cat(1, intDfBlock{:});
catch ME
end
close(h)
if ~isempty(ME)
    rethrow(ME);
end
intStats.chiSqSignificance = intChiSqSignificance;
intStats.chisq = intChiSq;
intStats.df = intDf;
chiSqSignificance = intChiSqSignificance(intIdVec,:);
chisq = intChiSq(intIdVec,:);
df= intDf(intIdVec,:);
stats.chiSqSignificance = chiSqSignificance;
stats.chisq = chisq;
stats.df = df;
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

