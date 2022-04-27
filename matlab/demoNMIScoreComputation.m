%% Calculate NMI score
addpath(genpath('AIDFUNCTIONS'));
DATASET='STAGE00DATA';
VAR_THRES=80;
LEVELS=4;
REDUCTION=1;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'LMU Serif');
cluster_path = @(component)sprintf('../results/%s/hierarchicalClusteringDemo/%s/%s_reduction%d/levels%d/segmentation.mat', ...
    component,DATASET, component, 1 / REDUCTION, LEVELS);
input_mat_path = @(component)sprintf('../results/%s/hierarchicalClusteringDemo/%s/%s_reduction%d/levels%d/input.mat', ...
    component, DATASET, component, 1 / REDUCTION, LEVELS);
asymmetry_cluster=load(cluster_path('asymmetry')).clusterArray;
symmetry_cluster=load(cluster_path('symmetry')).clusterArray;
load(input_mat_path('asymmetry'), 'preprocLandmarksIndices');
atlasName = 'Desikan_Killiany';
atlas = loadAtlas(atlasName,'L');
atlasIndices = repmat(atlas.index(preprocLandmarksIndices),1,5)';

OUTPUT_DIR = ['../results/asymmetryVsSymmetryNMI/' DATASET '/']; 
if ~isfolder(OUTPUT_DIR), mkdir(OUTPUT_DIR); end
%%
ret = cell(1,3);
for k=1:3
    switch k
        case 1
            cluster=asymmetry_cluster;
        case 2
            cluster=symmetry_cluster;
        case 3
            cluster = atlasIndices;
    end
    res = zeros(size(cluster,1)-1, size(cluster, 2));
    for row=2:size(cluster,1)
        clusterLevel = cluster(row,:);
        uniques = unique(clusterLevel);
        for u=1:length(uniques)
            flag = clusterLevel == uniques(u);
            res(row - 1, flag) = u;
        end
    end
    ret{k} = res;
end
%%
new_clustering1 = ret{1};
new_clustering2 = ret{2};
new_atlas = ret{3};
NMI_score = zeros(size(new_clustering1,1),3);
for level=1:size(new_clustering1,1)
    NMI_score(level, 1) = D_NMI(new_clustering1(level,:)', new_atlas(level,:)');
    NMI_score(level, 2) = D_NMI(new_clustering2(level,:)', new_atlas(level,:)');
    NMI_score(level, 3) = D_NMI(new_clustering1(level,:)', new_clustering2(level,:)');
end
%%
level = 1:size(new_clustering1,1);
level = level';
writetable(array2table(NMI_score, 'RowNames', string(level), 'VariableNames',  {'Asymmetry-Atlas','Symmetry-Atlas','Asymmetry-Symmetry'}), [OUTPUT_DIR 'result' num2str(1/REDUCTION)]);
%%
for level=1:4
    n = size(new_clustering1,2);
    n1 = length(unique(new_clustering1(level,:)));
    n2 = length(unique(new_clustering2(level,:)));
    n3 = length(unique(new_atlas(level,:)));
    
    mat = zeros(n1, ceil(n/n1));
    c = 1;
    r=1;
    f = ceil(n/n3);
    for i=1:n3
        if c + 2*f > size(mat,2)
            mat(r, c:c+f-1) = i;
            c=1;
            r = r+1;
        else
            mat(r, c: c + f-1) = i;
            c = c+f;
        end
    end
    c1 = reshape(mat', size(mat,1)*size(mat,2),1);
    c2 = reshape(repmat((1:size(mat,1)),size(mat,2),1),size(c1));
    
    max_score(level) = D_NMI(c1,c2);
end

%%
close(h)
%%

pNMI_score = NMI_score./max_score';

%%
level = 1:size(new_clustering1,1);
level = level';
writetable(array2table(pNMI_score, 'RowNames', string(level), 'VariableNames',  {'Asymmetry-Atlas','Symmetry-Atlas','Asymmetry-Symmetry'}), [OUTPUT_DIR 'result_scaled' num2str(1/REDUCTION)]);


