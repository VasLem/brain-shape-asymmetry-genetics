
function agg = aggregateExperimentIterations(nPicks,nSamplesPerPick, out,landmarksGroups, func)
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
                            for u=1:length(fs)
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
                            for u=1:length(fs)
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