function processDKAtlasAsymmetryAggregates(dkExpData_path, smallestPvalue)
    atlasName = 'Desikan_Killiany';
    side = 'L';
    [resultsDir,fName,~] = fileparts(dkExpData_path);
    data = load(dkExpData_path).data;
    x = cellfun(@struct2table, data.total,'uni',false);
    total = vertcat(x{:});
    total = total(2:end,:);
    
    atlas = loadAtlas(atlasName,side);
    indices = atlas.index;
    
    
    
    % data.values = oldData;
    availableToAtlasMask = indices ~= -1;
    
    c = categorical(indices(availableToAtlasMask));
    [~, gn] = grp2idx(c);
    atlasSideSpecifics = table(str2double(gn), 'VariableNames', {'ind'});
    
    
    
    for i=1:3
        for j=1:4
            colNames{i,j} = strcat(data.titleNames{i, 1}, '_', data.titleNames{i, j});
            if strcmp(data.titleNames{i, j}, 'Significance')
                atlasSideSpecifics.(colNames{i,j}) =  assignPvalueSignificance(total.(strcat("perm", data.titleNames{i,2})), smallestPvalue);
            elseif strcmp(data.titleNames{i, j}, 'p-Value')
                atlasSideSpecifics.(colNames{i,j}) =  total.(strcat("perm", data.titleNames{i,2}));
            else
                atlasSideSpecifics.(colNames{i,j}) =  total.(data.titleNames{i,j});
            end
        end
    end
    
    
    oriPart = table(indices(availableToAtlasMask), 'VariableNames',{'ind'});
    merged = join(oriPart,atlasSideSpecifics);
    for i=1:3
        for j=1:4
            ignoreValue = 0;
            if strcmp(data.titleNames{i, j}, 'p-value')
                ignoreValue = 1;
            end
            retAggregate = zeros(size(indices));
            retAggregate(:) = ignoreValue;
            retAggregate(availableToAtlasMask) = merged.(colNames{i,j});
            retData(i,j) = {retAggregate'};
        end
    end
    data.values = retData;
    %     save([resultsDir, 'atlas_', sides(s), atlasName, fil],'data','-v7');
    f = visualizeBrainAsymmetryData(data,[resultsDir 'aggregated_', fName]);
    set(f, 'Name', strcat([side, fName]));
    atlasSideSpecifics.label = atlas.labelsName(table2array(atlasSideSpecifics(:,1)) + 1)';
    resultForCSV = atlasSideSpecifics;
    writetable(resultForCSV, [resultsDir 'aggregated_', fName, '.csv']);
end




