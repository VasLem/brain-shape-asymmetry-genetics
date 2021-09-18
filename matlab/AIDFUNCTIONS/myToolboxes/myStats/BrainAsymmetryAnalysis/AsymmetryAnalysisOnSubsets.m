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
for s=1:length(nSamplesPerPick)
    samplePerPick=nSamplesPerPick(s);
    for iter=1:nPicks
        randIndex = randsample(n,samplePerPick);
        randX1 = X1(randIndex,:,:);
        randX2 = X2(randIndex,:,:);
        if ~isempty(landmarksGroups)
            uniqueGroups = unique(landmarksGroups);
            for g=1:length(uniqueGroups)
                mask = landmarksGroups == uniqueGroups(g);
                gMasks{g} = mask;
                repMask = repelem(mask, 3);
                gRet(g) = ProcrustesAnova2WayAsymmetryMEM(randX1(:, repMask,:),randX2(:, repMask,:),nIter,factor,nSplits);
            end
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
            clear gMasks gRet
            
        else
            out(s, iter) = ProcrustesAnova2WayAsymmetryMEM(randX1,randX2,nIter,factor,nSplits);
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
varargout{2} = aggregate(nPicks, nSamplesPerPick, out,landmarksGroups, 'mean');
if nOutputs == 2
    return
end
varargout{3} = aggregate(nPicks, nSamplesPerPick,out,landmarksGroups, 'std' );

end


function agg = aggregate(nPicks,nSamplesPerPick, out,landmarksGroups, func)
agg = struct();
fout = fieldnames(out(1,1));
for i = 1:length(fout)
    f = fout{i};
    if isstruct(out(1,1).(f))
        fds = fieldnames(out(1,1).(f));
        av = struct();
        for j=1:length(fds)
            attr = fds{j};
            for t = 1:length(nSamplesPerPick)
                for k=1:nPicks
                    v = out(t, k);
                    try
                        s(k, t,1)= v.(f).(attr);
                    catch
                        s(k, t,:)= v.(f).(attr);
                    end
                end
            end
            switch func
                case 'mean'
                    av.(attr) = squeeze(mean(s));
                case 'std'
                    av.(attr) = squeeze(std(s));
            end
        end
        agg.(f) = av;
    else
        for t = 1:length(nSamplesPerPick)
            for k=1:nPicks
                v = out(t,k);
                try
                    s(k,t,1)= v.(f);
                catch
                    try
                        s(k,t,:)= v.(f);
                    catch
                        s(k,t,:,:)= v.(f);
                    end
                end
            end
        end
        if isempty(landmarksGroups)
            switch func
                case 'mean'
                    agg.(f) = squeeze(mean(s));
                case 'std'
                    agg.(f) = squeeze(std(s));
            end
        else
            for i=1:length(unique(landmarksGroups))
                switch func
                    case 'mean'
                        try
                            agg.(f){i} = squeeze(mean(cell2mat(s(:,:,i)),1));
                        catch
                            p = cell2mat(s(:, :, i));
                            fs = fieldnames(p(1));
                            for u=1:length(f)
                                z = fs(u);
                                agg.(f){i}.(z{1}) = mean([p.(z{1})]);
                            end
                        end
                    case 'std'
                        try
                            agg.(f){i} = squeeze(std(cell2mat(s(:,:,i)),1));
                        catch
                            p = cell2mat(s(:, :, i));
                            fs = fieldnames(p(1));
                            for u=1:length(f)
                                z = fs(u);
                                agg.(f){i}.(z{1}) = std([p.(z{1})]);
                            end
                        end
                end
            end
        end
    end
    clear s;
end
end