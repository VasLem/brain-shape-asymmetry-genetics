function varargout  =  AsymmetryAnalysisOnSubsets(X1,X2,nSamplesPerPick, nPicks, nIter,factor,nSplits,landmarksGroups)
% nSamplesPerPick are the number of samples per subset
% nPicks are the number of picks
nOutputs = nargout;
varargout = cell(1,nOutputs);
n = size(X1, 1);
if nargin < 8, landmarksGroups = []; end
if nargin < 4, nPicks=1; end
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
    nGroups = 1;
end

[path,ID] = setupParForProgress(length(nSamplesPerPick) * nPicks * nGroups);
parfor s=1:length(nSamplesPerPick)
    samplePerPick=nSamplesPerPick(s);
    for iter=1:nPicks
        randIndex = randsample(n,samplePerPick);
        randX1 = X1(randIndex,:,:);
        randX2 = X2(randIndex,:,:);
        if ~isempty(landmarksGroups)
            gRet = cell(nGroups);
            
            for g=1:nGroups
                repMask = repelem(gMasks{g}, 3);
                gRet{g} = ProcrustesAnova2WayAsymmetryMEM(randX1(:, repMask,:),randX2(:, repMask,:),nIter,factor,nSplits, false);
                parfor_progress;
            end
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
            out(s, iter) = ProcrustesAnova2WayAsymmetryMEM(randX1,randX2,nIter,factor,nSplits, false);
            parfor_progress;
        end
       
    end
end
closeParForProgress(path,ID);
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

