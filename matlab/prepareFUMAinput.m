close all;clear;
restoredefaultpath;
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
DATA_DIR = '../SAMPLE_DATA/';


THREADS= 8;
MAX_NUM_FEATS = 0;
DATASET_INDEX = 1;

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
N_PARTITIONS = 31;
wholeTabs = cell(N_PARTITIONS,1);
for CHR=1:22
    disp(['CHR:' , num2str(CHR)]);
    CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
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


