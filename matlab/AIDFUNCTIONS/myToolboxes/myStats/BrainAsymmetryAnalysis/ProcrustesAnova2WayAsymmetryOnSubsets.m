function varargout  = ProcrustesAnova2WayAsymmetryOnSubsets(X1,X2,nSamplesPerPick, nPicks, t,factor,nSplits)
% nSamplesPerPick are the number of samples per subset
% nPicks are the number of picks
nOutputs = nargout;
varargout = cell(1,nOutputs);
n = size(X1, 1);

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
        
        out(s, iter) = ProcrustesAnova2WayAsymmetryMEM(randX1,randX2,t,factor,nSplits);
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
varargout{2} = aggregate(nPicks, nSamplesPerPick, out, 'mean');
if nOutputs == 2
    return
end
varargout{3} = aggregate(nPicks, nSamplesPerPick,out, 'std');

end


function agg = aggregate(nPicks,nSamplesPerPick, out,func)
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
        switch func
            case 'mean'
                agg.(f) = squeeze(mean(s));
            case 'std'
                agg.(f) = squeeze(std(s));
        end
    end
    clear s;
end
end