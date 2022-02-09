close all;clear;
restoredefaultpath;
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
DATA_DIR = '../SAMPLE_DATA/';


THREADS= 8;
MAX_NUM_FEATS = 0;
DATASET_INDEX = 2;

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
RESULTS_DIR = ['../results/genomeDemo/' DATASET_NAME '/'];
AVAILABLE_CHRS = [];
wholeTab = nan;
for CHR=1:22
    disp(['CHR:' , num2str(CHR)]);
    CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
    
    try
        tab = readtable([CHR_DIR 'chisq_stats.csv']);
        tab = tab(:,{'RSID', 'N', 'P_PAR1'});
        tab = renamevars(tab,{'RSID','P_PAR1'},{'rsID','P-value'});
        tab.chromosome = CHR * ones(height(tab),1);
        if ~isa(wholeTab, 'table')
            wholeTab = tab;
        else
            wholeTab = cat(1,wholeTab, tab);
        end
    catch
    end
end
writetable(wholeTab, [RESULTS_DIR 'noPartWholeCCA.csv']);