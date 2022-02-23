close all;clear;
restoredefaultpath;
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
DATA_DIR = '../SAMPLE_DATA/';

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
    inp_file = [CHR_DIR 'PartitionedGTLWoutBC_feats0significant_intervals.csv'];
    if ~isfile(inp_file)
        continue
    end
    ftab = readtable(inp_file);

    for partition=1:N_PARTITIONS
        tab = ftab(table2array(ftab(:,'PARTITION')) == partition,:);
        if isempty(tab)
            continue
        end
        tab = tab(:,{'CHR_ST', 'POS_ST', 'POS_EN'});
        tab.CHR_ST = cellstr(repmat(['chr' num2str(CHR)], height(tab), 1));
        fname =  [RESULTS_DIR sprintf('SigIntPart%02d.bed',partition)];
        if CHR ~= 1
            writetable(tab,fname,  "WriteMode","append",'FileType','text','WriteRowNames',false,'WriteVariableNames',false,'Delimiter',' ');
        else
            writetable(tab,fname,'FileType','text', 'WriteRowNames',false, 'WriteVariableNames',false, 'Delimiter',' ');
        end
    end
end


