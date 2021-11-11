function [stats, intStats] = runCCA(phenoPart, geno, intervals, intIdVec)
arguments
    phenoPart single
    geno
    intervals int32
    intIdVec int32
end
intChiSqSignificance = zeros(size(intervals, 1),1);

[...%~,~,
    ~,~, Q2, T22, perm2, rankY] = vl_mycanoncorr((1:size(phenoPart,1))', phenoPart);

if nargin < 4
    intIdVec = zeros(intervals(end,2),1);
    intIdVec(intervals(:,1)) = 1;
    intIdVec = cumsum(intIdVec);
end
if ~iscell(geno)
    geno = splitapply( @(x){x'}, geno', intIdVec );
end

[path,ID] = setupParForProgress(length(intervals));
parfor i=1:length(intervals)
%     tic;
    genoPart = double(geno{i});
    validSelection = ~any(isnan(genoPart), 2);
    [...%A,~,
        ~,st] = vl_mycanoncorr(genoPart(validSelection, :), phenoPart(validSelection, :), Q2(validSelection, :), T22, perm2, rankY);

    intChiSqSignificance(i) = st.pChisq(1);
%     intCoeffs{i} = A(:,1);
%     toc;
    parfor_progress;
end
closeParForProgress(path,ID);
intStats.chiSqSignificance = intChiSqSignificance;
% intStats.coeffs = intCoeffs;

chiSqSignificance = intChiSqSignificance(intIdVec);
% coeffs = intCoeffs{intIdVec};
stats.chiSqSignificance = chiSqSignificance;
% stats.coeffs = coeffs;
end