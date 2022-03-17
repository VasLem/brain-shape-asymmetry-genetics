close all;clear;
restoredefaultpath;
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));

% TO CHANGE
DATASET_INDEX = 1;
SUBSAMPLED_ID='not_subsampled';
IMPUTE_STRATEGY = 'beagle';

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
% CONSTANTS
DATA_DIR = '../SAMPLE_DATA/';
N_PARTITIONS = 31;
NO_PARTITION_THRES = 5*10^-8; % European in LD score

switch DATASET_INDEX
    case 1
        UKBIOBANK = 'UKBIOBANK';
        DATASET_NAME = 'STAGE00DATA';
        GENO_ID = 'sel19908';
    case 2
        UKBIOBANK = 'MY_UKBIOBANK';
        DATASET_NAME = 'BATCH2_2021_DATA';
        GENO_ID = 'sel16875_rmrel';
end
DATASET_ID = sprintf('%s/%s/%s', DATASET_NAME, IMPUTE_ID, SUBSAMPLED_ID);
LOADING_DIR = ['../results/genomeDemo/' DATASET_ID '/'];
RESULTS_DIR = ['../results/meta_analysis/' DATASET_ID '/'];
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
AVAILABLE_CHRS = [];
for CHR=1:22
    disp(['CHR:' , num2str(CHR)]);
    CHR_DIR = [LOADING_DIR 'chr' num2str(CHR) '/'];
    try
        ftab = readtable([CHR_DIR 'chisq_stats.csv']);
        for partition=1:N_PARTITIONS
            tab = ftab(:,{'RSID', 'N', 'A2', 'A1', ['P_PAR', num2str(partition)], ['CHI_PAR',num2str(partition)]});
            tab = renamevars(tab,{'RSID',['P_PAR' num2str(partition)], ['CHI_PAR', num2str(partition)]},{'rsID','P-value','ChiScore'});
            tab.chromosome = CHR * ones(height(tab),1);
            fname =  [RESULTS_DIR sprintf('CCAPart%02d.csv',partition)];
            if CHR ~= 1
                writetable(tab,fname,"WriteRowNames",CHR == 1, "WriteMode","append");
            else
                writetable(tab,fname);
            end
        end
    catch
    end
end
for partition=1:N_PARTITIONS
    fname =  [RESULTS_DIR sprintf('CCAPart%02d.csv',partition)];
    gzip(fname);
    delete(fname);
end


