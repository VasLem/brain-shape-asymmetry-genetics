function [stats, intStats] = runCCA(phenoPart, geno, intervals, intIdVec)
intChiSqSignificance = zeros(size(intervals, 1),1);
intChiSq = zeros(size(intervals, 1),1);
[~,~, Q2, T22, rankY] = vl_mycanoncorr((1:size(phenoPart,1))', phenoPart);
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

if ~exist('pagesvd','builtin') || ~exist('pagemtimes', 'builtin')
    [intChiSq,  intChiSqSignificance] = fallbackFunc(geno, phenoPart,Q2, T22, rankY);
else
    genoSizes = intervals(:,2) - intervals(:,1) + 1;
    flagWithNans = cellfun(@(x)any(any(isnan(x))), geno) == 1;
    
    for s=1:max(genoSizes)
    
    
        blockSize = 10000;
    
        flag = (genoSizes == s) & ~flagWithNans;
        total = sum(flag);
        part = geno(flag);
        blocksN = ceil(total/blockSize);
        f = waitbar(0,'1','Name',['Handling size ' num2str(s)]);
        clk = clock;
    
        intChiSqSignificanceSize = cell(blocksN,1);
        intChiSqSize =  cell(blocksN,1);
        for blockCnt=1: blocksN
            minInd = (1 + (blockCnt - 1) * blockSize);
            maxInd = min(total, (blockCnt * blockSize));
            block = part(minInd:maxInd);
            Xs = cat(3,block{:});            
            [intChiSqSize{blockCnt},  intChiSqSignificanceSize{blockCnt}] = vl_ccachisq1(Xs, phenoPart,Q2, T22, rankY);
            if blockCnt ==1
                is = etime(clock,clk);
                esttime = is * blocksN;
            end
            waitbar(blockCnt/blocksN,f,[num2str(blockCnt),'/',num2str(blocksN),', remaining time =',num2str(esttime-etime(clock,clk),'%4.1f'),'sec' ])
        end
        delete(f);
        intChiSqSignificance(flag) = cat(2, intChiSqSignificanceSize{:});
        intChiSq(flag) = cat(2, intChiSqSize{:});
    end
    
    % For cases with nans
    if any(flagWithNans)
        inds = find(flagWithNans);
        [intChiSq(inds), intChiSqSignificance(inds) ] = fallbackFunc(geno(inds));
        intChiSqSignificance(inds) = tmpSig;
    end
end

intStats.chiSqSignificance = intChiSqSignificance;
intStats.chisq = intChiSq;

chiSqSignificance = intChiSqSignificance(intIdVec);
chisq = intChiSq(intIdVec);
stats.chiSqSignificance = chiSqSignificance;
stats.chisq = chisq;
end


function [intChiSq, intChiSqSignificance] = fallbackFunc(geno, phenoPart, Q2, T22, rankY)
ppm = ParforProgressbar(length(geno), 'showWorkerProgress', true);
ME = [];
try
    intChiSq = ones(length(geno),1);
    intChiSqSignificance =  ones(length(geno),1);
    parfor i=1:length(geno)

        genoPart = single(geno);
        validSelection = ~any(isnan(genoPart), 2);
        [~,st] = vl_mycanoncorr(genoPart(validSelection, :), phenoPart(validSelection, :), Q2(validSelection, :), T22, rankY);
        intChiSq(i) =  st.pChisq(1);
        intChiSqSignificance(i) = st.chisq(1);
        ppm.increment();
    end
catch ME
end
delete(ppm);
if ~isempty(ME)
    rethrow(ME);
end

end
