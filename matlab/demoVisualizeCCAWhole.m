close all;clear;
restoredefaultpath;
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
DATA_DIR = '../SAMPLE_DATA/';


MODALITY = 'asymmetry'; %asymmetry, symmetry
MAX_NUM_FEATS = 0;
DATASET_INDEX = 0;
SUBSAMPLED = 0;
IMPUTE_STRATEGY = 'mean';
NO_PARTITION_THRES = 5*10^-8; % European in LD score

switch DATASET_INDEX
    case 0
        DATASET_NAME = 'joinedDatasets';
    case 1
        DATASET_NAME = 'STAGE00DATA';
    case 2
        DATASET_NAME = 'BATCH2_2021_DATA';
end
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

if SUBSAMPLED
    REDUCTION_ID = 'subsampled';
else
    REDUCTION_ID = 'not_subsampled';
end
RESULTS_ROOT = '../results/';
if DATASET_INDEX ~=0
    RESULTS_DIR = [RESULTS_ROOT, MODALITY '/genomeDemo/' DATASET_NAME '/' IMPUTE_ID '/' REDUCTION_ID '/'];
else
    RESULTS_DIR = [RESULTS_ROOT, MODALITY '/meta_analysis/' DATASET_NAME '/' IMPUTE_ID '/' REDUCTION_ID '/'];
   
end
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
%%
AVAILABLE_CHRS = [];
PNUM = 31;
part = cell(22, 1);
noPart = cell(22,1);
if DATASET_INDEX ~=0
    for CHR=1:22
        disp(['CHR:' , num2str(CHR)]);
        CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
        try
            load([CHR_DIR 'withPartCCA.mat'],'stats');
            part{CHR} = stats.chisqSignificance;
            noPart{CHR} = stats.chisqSignificance(:,1);
            AVAILABLE_CHRS(end+1) = CHR;
        catch
        end
    end
else
    partition = 1;
        PAR_FILE =  sprintf('%sCCAPart%02d.csv',RESULTS_DIR,partition);
       
        gunzip([PAR_FILE '.gz']);
        tab = readtable(PAR_FILE,"Delimiter",",","ReadVariableNames",1);
        
        delete(PAR_FILE);
        for CHR=1:22
                chrtab = tab(tab.chromosome == CHR,:);
                [~, a] = sort(chrtab.position);
                noPart{CHR} = table2array(chrtab(a, 'P_value')); 
                AVAILABLE_CHRS(end+1) = CHR;
        end 
end
disp(['Chromosomes correctly parsed:', num2str(AVAILABLE_CHRS)])
%%
f = figure;
f.Position = [100 100 1200 600];
hold on
cnt = 0;
xticks_pos = [];
st_edges = [];
en_edges = [];
for CHR=1:22
    if ~ismember(CHR, AVAILABLE_CHRS)
        continue
    end
    intStats = noPart{CHR};
    scatter(cnt + (1:length(intStats)), -log10(intStats),'.');
    xticks_pos(end + 1) = cnt + length(intStats) / 2;
    st_edges(end + 1) = cnt;
    cnt = cnt + length(intStats);
    en_edges(end + 1) = cnt;
end
xticks(xticks_pos)
xticklabels(arrayfun(@num2str, AVAILABLE_CHRS , 'UniformOutput', 0))
yline(-log10(NO_PARTITION_THRES));
ylabel('-log10p');
xlim([0, cnt]);
set(gca,'TickDir','out');
saveas(f, [RESULTS_DIR 'gwas.svg']);
close(f);
%%
if DATASET_INDEX ~= 0
pThres = NO_PARTITION_THRES;
bpThres = NO_PARTITION_THRES / 31;
fig = figure;
fig.Position = [100 100 1200 600];
hold on
cnt = 0;
xticks_pos = [];
st_edges = [];
en_edges = [];
for CHR=1:22
     if ~ismember(CHR, AVAILABLE_CHRS)
        continue
    end
    intStats = part{CHR};
    pNum = size(intStats, 2);
    for i=1:pNum
        signum = sum(intStats(:, i)<pThres);
        bsignum = sum(intStats(:, i)<bpThres);
        if signum > 0
            scatter(cnt + (1:length(intStats(:, i))), -log10(intStats(:, i)),'.', ...
                'DisplayName',['Chr. ' num2str(CHR) ', Part. ' num2str(i) ', # significant:', ...
                num2str(signum), ' (',num2str(bsignum),')']);
        end
    end
    xticks_pos(end + 1) = cnt + length(intStats(:, i)) / 2;
    st_edges(end + 1) = cnt;
    cnt = cnt + length(intStats(:, i));
    en_edges(end + 1) = cnt;
end

yline(-log10(pThres), '-k', 'DisplayName', 'Threshold');
yline(-log10(bpThres), '--k', 'DisplayName', '(Bonferroni Threshold)');


ylabel('-log10p');
xticks(xticks_pos)
xticklabels(arrayfun(@num2str, AVAILABLE_CHRS, 'UniformOutput', 0))
lgd = legend;
lgd.NumColumns = 4;
set(lgd, 'LimitMaxLegendEntries', false);
set(lgd,'Location','BestOutside');
ylabel('-log10p');
xlim([0, cnt]);
set(gca,'TickDir','out');
hold off;
saveas(fig, [RESULTS_DIR 'partitions_gwas.svg']);

close(fig);
end
