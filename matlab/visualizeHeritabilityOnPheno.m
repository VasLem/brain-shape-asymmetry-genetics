addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
DATASET_INDEX = 1;
switch DATASET_INDEX
    case 1
        UKBIOBANK = 'UKBIOBANK';
        DATASET = 'STAGE00DATA';
        GENO_ID = 'sel19908';
    case 2
        UKBIOBANK = 'MY_UKBIOBANK';
        DATASET = 'BATCH2_2021_DATA';
        GENO_ID = 'sel16875_rmrel';
end
LDSC_DIR = ['../results/ldsc/' DATASET '/'];
CLUSTER_DIR = ['../results/hierarchicalClusteringDemo/' DATASET '/'];
RESULTS_DIR = ['../results/visualizeHeritabilityOnPheno/' DATASET '/'];
clusterArray = load([CLUSTER_DIR  'asymmetry_reduction10/ccPriorSegmentation/levels4_mine/segmentation.mat']).clusterArray;
template = load([CLUSTER_DIR 'asymmetry_reduction10/ccPriorSegmentation/levels4_mine/input_info.mat']).preprocTemplate;
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
%%
heritabilility_file = [LDSC_DIR 'ldsc_heritability_results.csv'];
table = readtable(heritabilility_file,"ReadRowNames",   true,Delimiter='\t',ReadVariableNames=true);
%%
rows = table.Row;
for rowCnt=1:length(rows)
    row = rows{rowCnt};
    values = table2array(table(row,:));
    [fig, fig2, handles] = paintClusters(clusterArray, template, 4, false, values);
    close(fig)
    print(fig2, '-dpng', '-r300', [RESULTS_DIR, row, '.png'])
    height = fig2.Position(4);
    width = fig2.Position(3);
    close(fig2);
end
%%

in = cell(length(rowCnt),1);
for rowCnt = 1:length(rows)
    row = rows{rowCnt};
    in{rowCnt} = imread([RESULTS_DIR, row, '.png'],"png");
end
%%
finalFig = figure;
finalFig.Units='normalized';
finalFig.Position = [0.1,0.1, 4 * width, height];
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'none'); 
for rowCnt = 1:length(rows)
    row = rows{rowCnt};
    nexttile();
    imshow(in{rowCnt});
    title(row);
end
%%
saveas(finalFig,[RESULTS_DIR 'total.svg'])