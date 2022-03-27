close all; clear
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
RECOMPUTE_PARTS = 0;
SUBSAMPLED = 0;
MODALITY = 'symmetry';
DATASET_INDEX = 1;
IMPUTE_STRATEGY = 'mean';

switch DATASET_INDEX
    case 1
        DATASET = 'STAGE00DATA';
    case 2
        DATASET = 'BATCH2_2021_DATA';
end


switch IMPUTE_STRATEGY
    case 'no'
        IMPUTE_ID = 'not_imputed';
    case 'zero'
        IMPUTE_ID = 'zero_imputed';
    case 'median'
        IMPUTE_ID = 'median_imputed';
    case 'mean'
        IMPUTE_ID = 'mean_imputed';
    case 'beagle'
        IMPUTE_ID = 'beagle_imputed';
    otherwise
        error("IMPUTE_STRATEGY not understood, available options: no zero median beagle")
end

if SUBSAMPLED
    REDUCTION_ID='subsampled';
    REDUCTION=10;
else
    REDUCTION_ID='not_subsampled';
    REDUCTION=1;
end
RESULTS_ROOT = '../results/';
GENO_DIR = [RESULTS_ROOT, MODALITY '/genomeDemo/' DATASET '/' IMPUTE_ID '/' REDUCTION_ID '/'];

RESULTS_DIR = fullfile(pwd, [RESULTS_ROOT MODALITY '/visualizeCCAOnPheno/' DATASET '/' IMPUTE_ID '/' REDUCTION_ID '/']);

%%
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end

featMatsIds = ["WithBC", "WoutBC"];
featsClassesNames = cellstr(strcat('Chr', string(1:22)));
rescomputeParts = 0;


featMats = cell(2,1);
for i=1:2
    countsMat = zeros(31, 22);
    for j=1:31
        bCId = char(featMatsIds(i));
        for chr=1:22
            path = [ GENO_DIR 'chr' num2str(chr) '/PartitionedGTL' bCId 'significant_snps.csv'];
            if ~isfile(path)
                %                 disp(['Chromosome ' num2str(chr) ' snps file ' path ' not found. Skipping.']);
                continue
            end
            snpsTable = readtable(path);
            if isempty(snpsTable)
                continue
            end
            countsMat(j, chr) = sum((snpsTable.PARTITION == j) & (snpsTable.CHR==chr));
        end
    end
    countsMat(countsMat==0) = nan;
    featMats{i} = countsMat;
    
end
%%
drawFeaturesOnPolarPartitionsGraph(featMats, featMatsIds, featsClassesNames, MODALITY, RESULTS_DIR, REDUCTION, RECOMPUTE_PARTS)

