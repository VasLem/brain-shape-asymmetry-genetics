clear
close all
atlasName = 'Desikan_Killiany';
sides = ['L', 'R'];
experimentName = '19644_10_1000_3_1_100_100_3';
RESULTS_DIR = '../results/demo_asymmetry/';

mask = load('/usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat').index;
% oldData = load([RESULTS_DIR 'data_asymmetry.mat']).data;
% %%
% for i=1:3
%     for j=1:4
%         r = oldData(i,j);
%         r{1} = r{1}';
%         oldData(i,j) = r;
%     end
% end

for s=1:2
    atlas = load(['../SAMPLE_DATA/atlasses/', sides(s) , atlasName, '_Average.mat']).(atlasName);
    
    indices = atlas.index(mask);
    fil = ['data_' experimentName '.mat'];
    
    data = load([RESULTS_DIR fil]).data;
    % data.values = oldData;
    availableToAtlasMask = indices ~= -1;
    
    c = categorical(indices(availableToAtlasMask));
    [n, gn] = grp2idx(c);
    atlasSideSpecifics{s} = table(str2double(gn), 'VariableNames', {'ind'});
    

            
    for i=1:3
        for j=1:4
            k = 1;
            r = data.values(i,j);
            r = r{k}';
            colNames{i,j} = strcat(data.titleNames{i, 1}, '_', data.titleNames{i, j});
            func = @mean;
            if strcmp(data.titleNames{i, j}, 'p-value')
                func = @median;
            end
            atlasSideSpecifics{s}.(colNames{i,j}) = splitapply(func, r(availableToAtlasMask),n);
        end
    end
    oriPart = table(indices(availableToAtlasMask), 'VariableNames',{'ind'});
    merged = join(oriPart,atlasSideSpecifics{s});
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
%     data.values = retData;
%     save(['../results/demo_asymmetry/agg_', sides(s), atlasName, fil],'data','-v7');
%     f = visualizeBrainAsymmetryData(data,[RESULTS_DIR 'agg_', sides(s), atlasName, experimentName]);
%     set(f, 'Name', strcat([side(s), atlasName, experimentName]));
    atlasSideSpecifics{s}.label = atlas.labelsName(table2array(atlasSideSpecifics{s}(:,1)) + 1)';
end
resultForCSV =  [atlasSideSpecifics{1}; atlasSideSpecifics{2}];
writetable(resultForCSV, [RESULTS_DIR 'agg_', atlasName, experimentName, '.csv']);





