function varargout  =  AsymmetryAnalysisOnSubsets(X1,X2, nIter, nSamplesPerPick, nPicks, factor,nSplits,landmarksGroups, seed)
% X1: n_individuals x 3 x n_rep matrix
% X2: n_individuals x 3 x n_rep matrix
% nIter: the number of iterations to apply for the permutation test. If
% 0, no permutation test will be performed
% nSamplesPerPick: the number of samples per random subset picked
% nPicks: the number of random subsets picked
% factor: factor to divide X1 and X2 if they are integers before 2-way
% anova
% nSplits: number of splits for computational efficiency
% landmarksGroups: whether to group landmarks before applying 2-way anova
% seed: non negative integer for reproducibility of random number generation
arguments
X1
X2
nIter=0
nSamplesPerPick=size(X1,1)
nPicks=1
factor=10000
nSplits=1
landmarksGroups=[]
seed=-1
end
if seed>=0, rng(seed); end
nOutputs = nargout;
varargout = cell(1,nOutputs);
n = size(X1, 1);
if nargin < 3, nSamplesPerPick=n; end
if ~isvector(nSamplesPerPick)
    if (nSamplesPerPick == 0)
        nSamplesPerPick = n;
    end
    if (nSamplesPerPick > n)
        warning(["nSamplesPerPick(" num2str(nSamplesPerPick) ") greater than n(" num2str(n) "), changing it to be n"])
        nSamplesPerPick = n;
    end
    if ((nSamplesPerPick == n) && (nPicks >1))
        warning("nPicks reverted to 1, as nSamplesPerPick==n")
        nPicks = 1;
    end
end

if ~isempty(landmarksGroups)
    uniqueGroups = unique(landmarksGroups);
    nGroups = length(uniqueGroups);
    gMasks = cell(nGroups);
    for g=1:nGroups
        mask = landmarksGroups == uniqueGroups(g);
        gMasks{g} = mask;
        
    end
else
    gMasks = cell(1);
    uniqueGroups = 0; % req bc of parfor
    nGroups = 1; % req bc of parfor
end


for s=1:length(nSamplesPerPick)
    if length(nSamplesPerPick) > 1, disp(['Handling sample size ', num2str(s)]); end
    samplePerPick=nSamplesPerPick(s);
    for iter=1:nPicks
        disp(['Handling pick ',num2str(iter)])
        randIndex = randsample(n,samplePerPick);
        randX1 = X1(randIndex,:,:);
        randX2 = X2(randIndex,:,:);
        if ~isempty(landmarksGroups)
            gRet = cell(nGroups);
            ppb = ParforProgressbar(nGroups);
            for g=1:nGroups
                repMask = repelem(gMasks{g}, 3);
                gRet{g} = ProcrustesAnova2WayAsymmetryMEM(randX1(:, repMask,:),randX2(:, repMask,:),nIter,factor,nSplits, false);
                ppb.increment();
            end
            delete(ppb);
            gRet = [gRet{:}];
            fields = fieldnames(gRet(1).LM);
            gTret = struct;
            gTret.LM = struct;
            for f=1:length(fields)
                name = fields(f);
                name = name{1};
                gTret.LM.(name) = zeros(1, size(randX1,2)/3);
                for g=1:length(uniqueGroups)
                    t = gTret.LM.(name);
                    t(gMasks{g}) = gRet(g).LM.(name);
                    gTret.LM.(name) = t;
                end
            end
            meanX1 = cell(1,length(gRet));
            [meanX1{:}] = gRet.MeanX1;
            gTret.MeanX1 = meanX1;
            meanX2 = cell(1,length(gRet));
            [meanX2{:}] = gRet.MeanX2;
            gTret.MeanX2 = meanX2;
            meanDiff = cell(1,length(gRet));
            [meanDiff{:}] = gRet.MeanDiff;
            gTret.MeanDiff = meanDiff;
            total = cell(1,length(gRet));
            [total{:}] = gRet.Total;
            gTret.Total = total;
            out(s,iter) = gTret;
        else
            out(s, iter) = ProcrustesAnova2WayAsymmetryMEM(randX1,randX2,nIter,factor,nSplits, true);
           
        end
       
    end
end
varargout{1} = out;
if nOutputs == 1
    return
end

if nPicks == 1
    varargout{2} = out(1);
    return
end
varargout{2} = aggregateExperimentIterations(nPicks, nSamplesPerPick, out,landmarksGroups, 'mean');
if nOutputs == 2
    return
end
varargout{3} = aggregateExperimentIterations(nPicks, nSamplesPerPick,out,landmarksGroups, 'std' );

end

