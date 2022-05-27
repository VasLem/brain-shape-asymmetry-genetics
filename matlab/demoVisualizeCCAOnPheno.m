close all; clear
setuplatex
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
RECOMPUTE_PARTS = 1;
SUBSAMPLED = 0;
MODALITY = 'asymmetry';
DATASET_INDEX = 0;
IMPUTE_STRATEGY = 'mean';

switch DATASET_INDEX
    case 0
        DATASET = 'joinedDatasets';
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
if DATASET_INDEX ~= 0
    GENO_DIR = [RESULTS_ROOT, MODALITY '/genomeDemo/' DATASET '/' IMPUTE_ID '/' REDUCTION_ID '/'];
else
    GENO_DIR = [RESULTS_ROOT, MODALITY '/meta_analysis/' DATASET '/' IMPUTE_ID '/' REDUCTION_ID '/'];
end
RESULTS_DIR = fullfile(pwd, [RESULTS_ROOT MODALITY '/visualizeCCAOnPheno/' DATASET '/' IMPUTE_ID '/' REDUCTION_ID '/']);

%%
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end

featMatsIds = ["WithBC", "WoutBC"];
featsClassesNames = cellstr(strcat('Chr', string(1:22)));



featMats = cell(2,1);
for i=1:2
    countsMat = zeros(31, 22);
    parfor j=1:31
        bCId = char(featMatsIds(i));
        for chr=1:22
            if DATASET_INDEX ~= 0
                path = [ GENO_DIR 'chr' num2str(chr) '/PartitionedGTL' bCId 'significant_snps.csv'];
            else
                path = [ GENO_DIR '/' sprintf('CCAPart%02d.csv',j)];
                gunzip([path, '.gz']);
            end
            if ~isfile(path)
                %                 disp(['Chromosome ' num2str(chr) ' snps file ' path ' not found. Skipping.']);
                continue
            end
            if DATASET_INDEX ~= 0
                snpsTable = readtable(path);
            else
                snpsTable = readtable(path, "ReadVariableNames",1,Delimiter=',');
            end
            if isempty(snpsTable)
                continue
            end
            if DATASET_INDEX == 0
                if strcmp(bCId,'WoutBC')
                    t = 5e-8;
                else
                    t = 5e-8/31;
                end
                countsMat(j, chr) = sum((snpsTable.chromosome==chr) & snpsTable.P_value<=t);
            else
                countsMat(j, chr) = sum((snpsTable.PARTITION == j) & (snpsTable.CHR==chr));
            end
            if DATASET_INDEX == 0
                delete(path)
            end
        end
    end
    countsMat(countsMat==0) = nan;
    featMats{i} = countsMat;

end
%%

drawFeaturesOnPolarPartitionsGraph(featMats, featMatsIds, featsClassesNames, MODALITY, RESULTS_DIR, REDUCTION, RECOMPUTE_PARTS,  1)

