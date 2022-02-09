function [stats, intStats] = runCCA(phenoPart, geno, intervals, intIdVec)
intChiSqSignificance = zeros(size(intervals, 1),1);
intChiSq = zeros(size(intervals, 1),1);
% disp("Computing Q, R and the rank of the phenotype array..")
[...%~,~,
    ~,~, Q2, T22, rankY] = vl_mycanoncorr((1:size(phenoPart,1))', phenoPart);
% disp("Computing CCA..")
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
        bSize = maxInd - minInd + 1;
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
ppm = ParforProgressbar(length(inds), 'showWorkerProgress', true);
ME = [];
try
tmpSig = double(length(inds),1);
tmpSc =  double(length(inds),1);
  parfor i=1:length(inds)
 
        genoPart = single(geno{inds(i)});
        validSelection = ~any(isnan(genoPart), 2);
        [~,st] = vl_mycanoncorr(genoPart(validSelection, :), phenoPart(validSelection, :), Q2(validSelection, :), T22, rankY);
      tmpSig(i) =  st.pChisq(1);
      tmpSc(i) = st.chisq(1);
        ppm.increment();
   end
catch ME 
end
delete(ppm);
if ~isempty(ME)
   rethrow(ME);
end
intChiSqSignificance(inds) = tmpSig;
intChiSq(inds) = tmpSc;
end
intStats.chiSqSignificance = intChiSqSignificance;
intStats.chisq = intChiSq;
% intStats.coeffs = intCoeffs;

chiSqSignificance = intChiSqSignificance(intIdVec);
chisq = intChiSq(intIdVec);
% coeffs = intCoeffs{intIdVec};
stats.chiSqSignificance = chiSqSignificance;
stats.chisq = chisq;
% stats.coeffs = coeffs;
end