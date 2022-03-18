addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));

DATASET_INDEX = 1;
SUBSAMPLED=0;
IMPUTE_STRATEGY = 'mean';

switch IMPUTE_STRATEGY
    case 'no'
        IMPUTE_ID = 'not_imputed';
    case 'zero'
        IMPUTE_ID = 'zero_imputed';
    case 'median'
        IMPUTE_ID = 'median_imputed';
    case 'beagle'
        IMPUTE_ID = 'beagle_imputed';
    case 'mean'
        IMPUTE_ID = 'mean_imputed';
    otherwise
        error("IMPUTE_STRATEGY not understood, available options: no zero median beagle")
end


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
if SUBSAMPLED    
    SUBSAMPLED_ID = 'subsampled';
    REDUCTION=10;
else
    SUBSAMPLED_ID='not_subsampled';
    REDUCTION = 1;
end

DATASET_ID = sprintf('%s/%s/%s', DATASET_NAME, IMPUTE_ID, SUBSAMPLED_ID);

LDSC_DIR = ['../results/ldsc/' DATASET_ID '/h2/'];
CLUSTER_DIR = ['../results/hierarchicalClusteringDemo/' DATASET '/'];
RESULTS_DIR = ['../results/visualizeHeritabilityOnPheno/' DATASET_ID '/'];
load([CLUSTER_DIR  'asymmetry_reduction' num2str(REDUCTION) '/levels4/segmentation.mat'], 'clusterArray');
load([CLUSTER_DIR 'asymmetry_reduction' num2str(REDUCTION)  '/levels4/input_info.mat'], 'preprocTemplate');
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
%%
heritabilility_file = [LDSC_DIR 'h2_results.csv'];
table = readtable(heritabilility_file,"ReadRowNames",   true,Delimiter='\t',ReadVariableNames=true);
%%
rows = table.Row;
for rowCnt=1:length(rows)
    row = rows{rowCnt};
    values = table2array(table(row,:));
    [fig, fig2, handles] = paintClusters(clusterArray, preprocTemplate, 4, false, values);
    close(fig)
    saveas(fig2, [RESULTS_DIR, row, '.svg'])
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