DATASET_INDEX = 1;

IMPUTE_ID = 'median_imputed';
SUBSAMPLED = 0;

if SUBSAMPLED
    REDUCTION_ID='subsampled';
else
    REDUCTION_ID='not_subsampled';
end

switch DATASET_INDEX
    case 1
        DATASET_NAME = 'STAGE00DATA';
    case 2
        DATASET_NAME = 'BATCH2_2021_DATA';
end
RESULTS_ROOT = '../results/';

INPUT_DIR = [RESULTS_ROOT, 'genomeDemo/' DATASET_NAME '/' IMPUTE_ID '/' REDUCTION_ID '/'];
OUT_PATH = [RESULTS_ROOT, 'genomeDemo/' DATASET_NAME '/' IMPUTE_ID '/' REDUCTION_ID '/median_values.csv'];
whole_table = [];
for CHR=1:22
    load(populateCCAWorkspace(INPUT_DIR, INPUT_DIR, CHR));
    tab = readtable(IMPUTE_INFO_OUT,'Delimiter',' ');
    tab.CHR = zeros(height(tab),1) + CHR;
    if isempty(whole_table)
        whole_table = tab;
    else
        whole_table = [whole_table; tab];
    end
end
writetable(whole_table, OUT_PATH, 'Delimiter', ' ');