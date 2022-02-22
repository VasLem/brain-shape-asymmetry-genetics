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
covGenoPath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
disp("Initializing genomic analysis..")
try
    parpool('local',THREADS);
catch
end
%%
AVAILABLE_CHRS = [];
PNUM = 31;
part = cell(2, 1);

for CHR=1:22
    disp(['CHR:' , num2str(CHR)]);
    CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
    try
        load([CHR_DIR 'withPartCCA.mat'],'gTLPartStats');
        part{CHR} = gTLPartStats.chiSqSignificance;
        noPart{CHR} = gTLPartStats.chiSqSignificance(:,1);
        AVAILABLE_CHRS(end+1) = CHR;
    catch
    end
end
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
saveas(f, [RESULTS_DIR 'noPartWholeCCA_logPlot.svg']);
close(f);
%%
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
saveas(fig, [RESULTS_DIR 'gTLPartWholeCCA_logPlot.svg']);

close(fig);

