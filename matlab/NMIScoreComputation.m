%% Calculate NMI score
DATASET='STAGE00DATA';
VAR_THRES=80;
LEVELS=4;
REDUCTION=1;
cluster_path = @(component)sprintf('../results/hierarchicalClusteringDemo/%s/%s_reduction%d/levels%d/segmentation.mat', ...
    DATASET,component,1/REDUCTION,LEVELS);
asymmetry_cluster=load(cluster_path('asymmetry')).clusterArray;
symmetry_cluster=load(cluster_path('symmetry')).clusterArray;
%%
ret = cell(1,2);
for k=1:2
    if k ==1
        cluster=asymmetry_cluster;
    else
        cluster=symmetry_cluster;
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
NMI_score = zeros(size(new_clustering1,1),1);
for level=1:size(new_clustering1,1)
    NMI_score(level) = D_NMI(new_clustering1(level,:)', new_clustering2(level,:)');
end
level = 1:size(new_clustering1,1);
level = level';
writetable(table(NMI_score, level), ['../results/hierarchicalClusteringDemo/' DATASET '/asymmetryVsSymmetryNMI' num2str(1/REDUCTION)]);