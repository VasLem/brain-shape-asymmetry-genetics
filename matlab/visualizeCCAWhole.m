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

for CHR=1:22
    disp(['CHR:' , num2str(CHR)]);
    CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
    noPart{CHR} = load([CHR_DIR 'noPartCCA.mat']);
    part{CHR} = load([CHR_DIR 'withPartCCA.mat']);
end
%%
f = figure;
f.Position = [100 100 1200 600];
hold on
cnt = 0;

for CHR=1:22
    intStats = noPart{CHR}.noPartitionIntStats;
    scatter(cnt + (1:length(intStats.chiSqSignificance)), -log10(intStats.chiSqSignificance),'.');
    xticks_pos(CHR) = cnt + length(intStats.chiSqSignificance) / 2;
    st_edges(CHR) = cnt;
    cnt = cnt + length(intStats.chiSqSignificance);
    en_edges(CHR) = cnt;
end
xticks(xticks_pos)
xticklabels(arrayfun(@num2str, 1:22 , 'UniformOutput', 0))
yline(-log10(NO_PARTITION_THRES));
ylabel('-log10p');
xlim([0, cnt]);
set(gca,'TickDir','out');
saveas(f, [RESULTS_DIR 'noPartWholeCCA_logPlot.png']);
close(f);
%%
pThres = NO_PARTITION_THRES;
fig = figure;
fig.Position = [100 100 1200 600];
hold on
cnt = 0;
for CHR=1:22
    intStats = part{CHR}.gTLPartStats;
    pNum = size(intStats.chiSqSignificance, 1);
    for i=1:pNum
        signum = sum(intStats.chiSqSignificance(i, :)<pThres);
        if signum > 0
            scatter(cnt + (1:length(intStats.chiSqSignificance(i, :))), -log10(intStats.chiSqSignificance(i, :)),'.', ...
                'DisplayName',['Chr. ' num2str(CHR) ', Part. ' num2str(i) ', # significant:', ...
                num2str(signum)]);
        end
    end
    xticks_pos(CHR) = cnt + length(intStats.chiSqSignificance(i, :)) / 2;
    st_edges(CHR) = cnt;
    cnt = cnt + length(intStats.chiSqSignificance(i, :));
    en_edges(CHR) = cnt;
end

yline(-log10(pThres), 'DisplayName', 'Threshold');

ylabel('-log10p');
xticks(xticks_pos)
xticklabels(arrayfun(@num2str, 1:22 , 'UniformOutput', 0))
lgd = legend;
lgd.NumColumns = 4;
set(lgd, 'LimitMaxLegendEntries', false);
set(lgd,'Location','BestOutside');
ylabel('-log10p');
xlim([0, cnt]);
set(gca,'TickDir','out');
hold off;
saveas(fig, [RESULTS_DIR 'gTLPartWholeCCA_logPlot.png']);
close(fig);

%%
pThres = NO_PARTITION_THRES / 31;
fig = figure;
fig.Position = [100 100 1200 600];
hold on
cnt = 0;
for CHR=1:22
    intStats = part{CHR}.gTLPartStats;
    pNum = size(intStats.chiSqSignificance, 1);
    for i=1:pNum
        signum = sum(intStats.chiSqSignificance(i, :)<pThres);
        if signum > 0
            scatter(cnt + (1:length(intStats.chiSqSignificance(i, :))), -log10(intStats.chiSqSignificance(i, :)),'.', ...
                'DisplayName',['Chr. ' num2str(CHR) ', Part. ' num2str(i) ', # significant:', ...
                num2str(signum)]);
        end
    end
    xticks_pos(CHR) = cnt + length(intStats.chiSqSignificance(i, :)) / 2;
    st_edges(CHR) = cnt;
    cnt = cnt + length(intStats.chiSqSignificance(i, :));
    en_edges(CHR) = cnt;
end

yline(-log10(pThres), 'DisplayName', 'Threshold');

ylabel('-log10p');
xticks(xticks_pos)
xticklabels(arrayfun(@num2str, 1:22 , 'UniformOutput', 0))
lgd = legend;
lgd.NumColumns = 4;
set(lgd, 'LimitMaxLegendEntries', false);
set(lgd,'Location','BestOutside');
ylabel('-log10p');
xlim([0, cnt]);
set(gca,'TickDir','out');
hold off;
saveas(fig, [RESULTS_DIR 'gTLPartWholeCCABonf_logPlot.png']);
close(fig);




