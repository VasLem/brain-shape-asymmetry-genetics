addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
setuplatex;
DATASET_INDEX = 0;
SUBSAMPLED=0;
IMPUTE_STRATEGY = 'mean';
MODALITY = 'asymmetry'; %asymmetry,symmetry

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
    case 0
        DATASET = 'joinedDatasets';
    case 1
        DATASET = 'STAGE00DATA';
    case 2
        DATASET = 'BATCH2_2021_DATA';
end
if SUBSAMPLED    
    SUBSAMPLED_ID = 'subsampled';
    REDUCTION=10;
else
    SUBSAMPLED_ID='not_subsampled';
    REDUCTION = 1;
end

DATASET_ID = sprintf('%s/%s/%s', DATASET, IMPUTE_ID, SUBSAMPLED_ID);

LDSC_DIR = ['../results/' MODALITY '/ldsc/' DATASET_ID '/h2/'];
CLUSTER_DIR = ['../results/' MODALITY '/hierarchicalClusteringDemo/STAGE00DATA/'];
RESULTS_DIR = ['../results/' MODALITY '/visualizeHeritabilityOnPheno/' DATASET_ID '/'];
load([CLUSTER_DIR  MODALITY '_reduction' num2str(REDUCTION) '/levels4/segmentation.mat'], 'clusterArray');
load([CLUSTER_DIR MODALITY '_reduction' num2str(REDUCTION)  '/levels4/input_info.mat'], 'preprocTemplate');
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
%%
heritabilility_file = [LDSC_DIR 'h2_results.csv'];
table = readtable(heritabilility_file,"ReadRowNames",   true,Delimiter='\t',ReadVariableNames=true);
%%
rows = table.Row;
for rowCnt=1:length(rows)
    row = rows{rowCnt};
    values = table2array(table(row,:));
    map = repmat(linspace(0.2,1,100), [3,1])';

    [fig, fig2, handles] = paintClusters(clusterArray, preprocTemplate, 4, false, values,'k',map);
    close(fig)
    set(fig2, 'InvertHardCopy', 'off');
    set(fig2, 'Color', 'white');

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
    if rowCnt == 3
        row = ['$\mathrm{' row '}$'];
    end
    title(row);
end
%%
saveas(finalFig,[RESULTS_DIR 'total.svg'])