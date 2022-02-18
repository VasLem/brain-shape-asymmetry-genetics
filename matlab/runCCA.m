function [stats, intStats] = runCCA(pheno, geno, intervals, intIdVec)
blockSize = 10000;
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
genoBlocks = cell(blocksN,1);
for blockCnt=1: blocksN
    minInd = (1 + (blockCnt - 1) * blockSize);
    maxInd = min(total, (blockCnt * blockSize));
    genoBlocks{blockCnt} = geno(minInd:maxInd);
end

ppm = ParforProgressbar(blocksN, 'showWorkerProgress', true);
ME = [];
try
    parfor blockCnt=1: blocksN
        minInd = (1 + (blockCnt - 1) * blockSize);
        maxInd = min(total, (blockCnt * blockSize));
        currBlockSize = maxInd - minInd + 1;
        genoBlock = genoBlocks{blockCnt};
        finite_mask = cellfun(@(x)all(x~=255,2), genoBlock,'UniformOutput',false);
        finite_mask = cat(2,finite_mask{:})';
        [patterns, ~, patternsInds] = unique(finite_mask,'rows');
        fallback = ~exist('pagesvd','builtin') || ~exist('pagemtimes', 'builtin');
        blockSizes = genoSizes(minInd:maxInd);
        intChiSqSignificancePattern = ones(currBlockSize, phenoPartNum);
        intChiSqPattern = zeros(currBlockSize,phenoPartNum);
        intDfPattern = zeros(currBlockSize,phenoPartNum);
        for pcnt=1:size(patterns,1)
            pattern = patterns(pcnt,:);
            patternFlag= patternsInds == pcnt;
            patternSizes = blockSizes(patternFlag);
            patternNSnps = sum(patternFlag);
            for phenoPartInd=1:length(pheno)
                phenoPart = pheno{phenoPartInd};
                phenoPart = phenoPart(pattern, :);

                genoPattern =cellfun(@(x)x(pattern), genoBlock(patternFlag), 'UniformOutput', false);
                if patternNSnps==1
                        st = vl_mycanoncorr(double(genoPattern{1}), phenoPart);
                        intChiSqPattern(patternFlag, phenoPartInd)=  st.chisq;
                        intChiSqSignificancePattern(patternFlag, phenoPartInd) = st.pChisq;
                        intDfPattern(patternFlag, phenoPartInd) = st.df1;
                    continue
                end
                
                [~, Q2, T22, rankY] = vl_mycanoncorr([], phenoPart);
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
        ppm.increment();
    end
    intChiSqSignificance = cat(1, intChiSqSignificanceBlock{:});
    intChiSq = cat(1, intChiSqBlock{:});
    intDf = cat(1, intDfBlock{:});
catch ME
end
delete(ppm);
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
