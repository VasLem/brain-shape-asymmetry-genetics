clear
close all
atlasName = 'Desikan_Killiany';
side = 'L';
experimentName = '19644_10_1000_3_1_100_100_3_DK';
RESULTS_DIR = '../results/demo_asymmetry/';
addpath(genpath('AIDFUNCTIONS'));

% oldData = load([RESULTS_DIR 'data_asymmetry.mat']).data;
% %%
% for i=1:3
%     for j=1:4
%         r = oldData(i,j);
%         r{1} = r{1}';
%         oldData(i,j) = r;
%     end
% end
fil = ['data_' experimentName '.mat'];
data = load([RESULTS_DIR fil]).data;
x = cellfun(@struct2table, data.total,'uni',false);
total = vertcat(x{:});
total = total(2:end,:);

atlas = loadAtlas(atlasName,side);
indices = atlas.index;



% data.values = oldData;
availableToAtlasMask = indices ~= -1;

c = categorical(indices(availableToAtlasMask));
[n, gn] = grp2idx(c);
atlasSideSpecifics = table(str2double(gn), 'VariableNames', {'ind'});



for i=1:3
    for j=1:4
        k = 1;
        r = data.values(i,j);
        r = r{k}';
        colNames{i,j} = strcat(data.titleNames{i, 1}, '_', data.titleNames{i, j});
        if strcmp(data.titleNames{i, j}, 'Significance')
            atlasSideSpecifics.(colNames{i,j}) =  assignPvalueSignificance(total.(strcat("perm", data.titleNames{i,2})));
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
%     save(['../results/demo_asymmetry/agg_', sides(s), atlasName, fil],'data','-v7');
f = visualizeBrainAsymmetryData(data,[RESULTS_DIR 'agg_', side, atlasName, experimentName]);
set(f, 'Name', strcat([side, atlasName, experimentName]));
atlasSideSpecifics.label = atlas.labelsName(table2array(atlasSideSpecifics(:,1)) + 1)';

resultForCSV = atlasSideSpecifics;
writetable(resultForCSV, [RESULTS_DIR 'agg_', atlasName, experimentName, '.csv']);





