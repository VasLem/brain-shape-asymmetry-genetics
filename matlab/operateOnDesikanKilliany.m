clear
close all
atlas = 'Desikan_Killiany';
lDesikan = load(['../SAMPLE_DATA/atlasses/L' atlas '_Average.mat']).(atlas);
rDesikan = load(['../SAMPLE_DATA/atlasses/R' atlas '_Average.mat']).(atlas);
mask = load('/usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat').index;
lIndices = lDesikan.index(mask);
rIndices = rDesikan.index(mask);
experimentName = '1965_10_1000_5_1_200_200_3';
fil = ['data_' experimentName '.mat'];
RESULTS_DIR = '../results/demo_asymmetry/';
data = load([RESULTS_DIR fil]).data;
indices = rIndices;
availableToAtlasMask = indices ~= -1;
%%
oldData = load([RESULTS_DIR 'data_asymmetry.mat']).data;
%%
for i=1:3
    for j=1:4
        r = oldData(i,j);
        r{1} = r{1}';
        oldData(i,j) = r;
    end
end
data.values = oldData;
%%images
for i=1:3
    for j=1:4
        k = 1;
        ignoreValue = 0;
        if strcmp(data.titleNames{i, j}, 'p-value')
            ignoreValue = 1;
        end
        r = data.values(i,j);
        r = r{k}';
        c = categorical(indices(availableToAtlasMask));
        [n, gn] = grp2idx(c);
        aggCol = splitapply(@mean, r(availableToAtlasMask),n);
        oriPart = table(r(availableToAtlasMask), indices(availableToAtlasMask), 'VariableNames',{'value','ind'});
        retPart = table(aggCol, str2double(gn), 'VariableNames',{'aggregate','ind'});
        merged = join(oriPart,retPart);
        retAggregate = zeros(size(indices));
        retAggregate(:) = ignoreValue;
        retAggregate(availableToAtlasMask) = merged.('aggregate');
        retData(i,j) = {retAggregate'};
    end
end
%%
data.values = retData;
%%
save(['../results/demo_asymmetry/agg_' fil],'data','-v7');

%%

f = visualizeBrainAsymmetryData(data,[RESULTS_DIR 'agg_results_asymmetry_r']);
