%ENV TO SET:
%%%%%%%%%%%%
%DATA_ROOT: the directory of the DATA (../SAMPLE_DATA/)
%DATASET_INDEX: dataset to use, 1 or 2 (1)
%RESULTS_ROOT: the directory of the results (../results/)
%SCRATCH_ROOT: the directory  of the temporary results (../results/)
%THREADS: Number of threads to use (max set by local)
%CHROMOSOME: Chomosome to analyze (1:22)
%BLOCK_SIZE: Block Size for CCA (2000)
%IMPUTE_STRATEGY: Whether to perform imputation (zero,median or beagle) or not (no). No imputation is quite slow. (no)
%SUBSAMPLED: Whether to use the subsampled phenotype, if no PHENO_PATH is
%provided. It modifies the saving and scratch directories (defaults to 0)
%PHENO_PATH: Whether to use a specific path for the mat file of the phenotype (defaults to the path where hierarchical clustering algorithm places it)
%%%%%%%%%%%%
close all;clear;
if ~isdeployed
    restoredefaultpath;
    addpath(genpath('.'));
    addpath(genpath('AIDFUNCTIONS'));
    addpath(genpath('../SNPLIB/'));
end
%%
UPDATE_FIGS = 1;
NO_PARTITION_THRES = 5*10^-8; % European in LD score
DEFAULT_CHRS = 1:22;
DEFAULT_SUBSAMPLED = 0;
DEFAULT_DATASET_INDEX = 1;
DEFAULT_MEDIAN_IMPUTE = 'zero';
DEFAULT_BLOCK_SIZE =2000;

SUBSAMPLED = getenv('SUBSAMPLED');
if(isempty(SUBSAMPLED))
    SUBSAMPLED = DEFAULT_SUBSAMPLED;
end
if SUBSAMPLED
    REDUCTION_ID='subsampled';
else
    REDUCTION_ID='not_subsampled';
end

IMPUTE_STRATEGY = getenv('IMPUTE_STRATEGY');
if(isempty(IMPUTE_STRATEGY))
    IMPUTE_STRATEGY = DEFAULT_MEDIAN_IMPUTE;
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
    otherwise
        error("IMPUTE_STRATEGY not understood, available options: no zero median beagle")
end

BLOCK_SIZE = getenv('BLOCK_SIZE');
if(isempty(BLOCK_SIZE))
    BLOCK_SIZE = DEFAULT_BLOCK_SIZE;
end


DATA_DIR = getenv('DATA_ROOT');
if(isempty(DATA_DIR))
    DATA_DIR = '../SAMPLE_DATA/';
end

SCRATCH_ROOT = getenv('SCRATCH_ROOT');
if(isempty(SCRATCH_ROOT))
    SCRATCH_ROOT = '../results/';
end


THREADS = getenv('THREADS');
if(isempty(THREADS))
    THREADS = parcluster('local');
    THREADS = THREADS.NumWorkers;
else
    if ~isnumeric(THREADS)
        THREADS=str2double(THREADS);
    end
end

CHRS = getenv("CHROMOSOME");
if(isempty(CHRS))
    CHRS = DEFAULT_CHRS;
else
    if ~isnumeric(CHRS)
        CHRS=str2double( strsplit(CHRS,','));
    end
end


DATASET_INDEX = getenv("DATASET_INDEX");
if (isempty(DATASET_INDEX))
    DATASET_INDEX = DEFAULT_DATASET_INDEX;
else
    disp(DATASET_INDEX)
    if ~isnumeric(DATASET_INDEX)
        DATASET_INDEX = str2double(DATASET_INDEX);
    end
end

disp(['Number of threads: ', num2str(THREADS)])
disp(['Location of data: ', DATA_DIR])
disp(['Location of temporary data:', SCRATCH_ROOT])
disp(['Median Imputation: ', num2str(IMPUTE_STRATEGY)])
disp(['Block Size: ', num2str(BLOCK_SIZE)]);
disp(['Using dataset: ', num2str(DATASET_INDEX)])
disp(['Chromosome analyzed: ', num2str(CHRS)])

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

RESULTS_ROOT = getenv('RESULTS_ROOT');
if(isempty(RESULTS_ROOT))
    RESULTS_ROOT = '../results/';
end
RESULTS_DIR = [RESULTS_ROOT, 'genomeDemo/' DATASET_NAME '/' IMPUTE_ID '/' REDUCTION_ID '/'];
SCRATCH_DIR = [SCRATCH_ROOT, 'genomeDemo/' DATASET_NAME '/' IMPUTE_ID '/' REDUCTION_ID '/'];
disp(['Location of results: ', RESULTS_DIR]);
disp(['Location of temporary results:', SCRATCH_DIR])
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
if ~isfolder(SCRATCH_DIR), mkdir(SCRATCH_DIR); end

COV_GENO_PATH = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];
PHENO_PATH = getenv('PHENO_PATH');
if(isempty(PHENO_PATH))
    PHENO_PATH = [RESULTS_ROOT, 'hierarchicalClusteringDemo/' DATASET_NAME '/asymmetry_reduction1/levels4/phenotype_varThres80.mat'];
end

disp(['Loading phenotype from: ', PHENO_PATH])

load(PHENO_PATH, 'PHENO','PHENO_IID');

%%
disp("Initializing genomic analysis..")
try
    parpool(THREADS);
catch
end

for CHR_IND=1:length(CHRS)
    tic;
    clearvars -except -regexp ^[A-Z|_]+$
    CHR = CHRS(CHR_IND);
    GENE_SET_METHOD = 'perSNP';
    load(populateCCAWorkspace(RESULTS_DIR, SCRATCH_DIR, CHR));
    if WITH_PART_CCA_PROC

        display(['Preprocessing chromosome ', num2str(CHR)])
        if strcmp(IMPUTE_STRATEGY,'beagle')
            GENO_PATH = [DATA_DIR 'IMAGEN/BRAIN/IMPUTED_GENOTYPES/' DATASET_NAME '/ukb_img_maf0.01_geno0.5_hwe1e-6_' GENO_ID '_beagle_chr' num2str(CHR)];
        else
            GENO_PATH = [DATA_DIR 'IMAGEN/BRAIN/' UKBIOBANK '/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_' GENO_ID '_chr' num2str(CHR)];
        end
        preprocessGenome(PHENO_IID, CHR, GENO_PATH, IMPUTE_STRATEGY, BLOCK_SIZE, RESULTS_DIR, SCRATCH_GENE_DIR, THREADS, SCRATCH_DIR);
        clear genoPruned
        toc;
    end
end
%%
for CHR_IND=1:length(CHRS)
    clearvars -except -regexp ^[A-Z|_]+$
    CHR = CHRS(CHR_IND);
    load(populateCCAWorkspace(RESULTS_DIR, SCRATCH_DIR, CHR));
    if WITH_PART_CCA_PROC
        disp(['Running CCA on chromosome ' , num2str(CHR)]);
        tic;
        fileDir = SCRATCH_GENE_DIR;
        outFiles = runCCA(PHENO,  SCRATCH_GENE_DIR, SCRATCH_OUTPUT_GENE_DIR, ~strcmp(IMPUTE_STRATEGY, 'no'));
        rmdir(SCRATCH_GENE_DIR);
        toc;
    end
end
%%
for CHR_IND=1:length(CHRS)
        clearvars -except -regexp ^[A-Z|_]+$
    CHR = CHRS(CHR_IND);
    load(populateCCAWorkspace(RESULTS_DIR, SCRATCH_DIR, CHR));
    if WITH_PART_CCA_PROC
        disp(['Collecting results of chromosome ' , num2str(CHR)]);
        tic;
        blocksN = length(dir(SCRATCH_OUTPUT_GENE_DIR)) - 2; % removing . and ..
        intchisqSignificanceBlock = cell(blocksN,1);
        intChiSqBlock =  cell(blocksN,1);
        intDfBlock =  cell(blocksN,1);
        for blockCnt=1:blocksN
            load([SCRATCH_OUTPUT_GENE_DIR, num2str(blockCnt), '.mat']);
            intchisqSignificanceBlock{blockCnt} = chisqSignificance;
            intChiSqBlock{blockCnt} = chisq;
            intDfBlock{blockCnt} = df;
        end
        intStats.chisqSignificance = cat(1, intchisqSignificanceBlock{:});
        intStats.chisq = cat(1, intChiSqBlock{:});
        intStats.df = cat(1, intDfBlock{:});
        load(PLINK_DATA_INFO_OUT, 'intervals');
        intIdVec = zeros(intervals(end,2),1);
        intIdVec(intervals(:,1)) = 1;
        intIdVec = cumsum(intIdVec);
        stats.chisqSignificance = intStats.chisqSignificance(intIdVec,:);
        stats.chisq = intStats.chisq(intIdVec,:);
        stats.df= intStats.df(intIdVec,:);
        toc;
        disp("Saving CCA results..")
        tic;
        save(WITH_PART_CCA_OUT, 'stats', 'intStats', '-v7.3');
        toc;
    end
end
%%
for CHR_IND=1:length(CHRS)
    clearvars -except -regexp ^[A-Z|_]+$
    CHR = CHRS(CHR_IND);
    load(populateCCAWorkspace(RESULTS_DIR, SCRATCH_DIR, CHR));
    disp(['Processing CCA results of chromosome ' , num2str(CHR)]);
    disp("Loading computed variables..")
    tic;
    load(WITH_PART_CCA_OUT,  'stats', 'intStats');
    load(PLINK_DATA_INFO_OUT, 'snpsPruned', 'intervals', 'sampleSizes');
    toc;
    BONFERRONI_THRES = NO_PARTITION_THRES / length(PHENO);
    if UPDATE_FIGS
        disp("Plotting results..")
        tic;
        plotSimpleGWAS(intervals, intStats.chisqSignificance(:, 1), CHR, NO_PARTITION_THRES,  [CHR_DIR  'noPartition']);
        plotPartitionsGWAS(intervals, intStats, CHR, BONFERRONI_THRES, NO_PARTITION_THRES, [CHR_DIR 'PartitionedGTL']);
        toc;
    end
    %% Significant SNPs tables extraction
    disp("Loading snpsPruned..")
    disp("Extracting significant SNPs tables..")
    tic;
    prepareSignificantTablesOnEachPartition(snpsPruned,  intStats , intervals, BONFERRONI_THRES, [CHR_DIR, 'PartitionedGTLWithBC']);
    prepareSignificantTablesOnEachPartition(snpsPruned,  intStats , intervals, NO_PARTITION_THRES, [CHR_DIR, 'PartitionedGTLWoutBC']);
    toc;
    %% LD regression csv
    disp("Extracting partition specific signigicant SNPs tables..")
    tic;
    saveLDRegressionTablesOnEachPartition(snpsPruned, PHENO, sampleSizes, intStats.chisq, intStats.chisqSignificance, intervals, CHR_DIR);
    tic;
    %
    disp("End of computation.")
    %%
    close all
    toc;
end

function output = saveLDRegressionTablesOnEachPartition(snpsPruned,  PHENO, sample_sizes, partScores, partSignificances, intervals, save_dir)
pNum = size(partSignificances ,2);
names = snpsPruned.Properties.VariableNames;
new_names = strrep(names,'ALT', 'A2');
new_names = strrep(new_names,'REF', 'A1');
output = snpsPruned;
output.Properties.VariableNames = new_names;
output = output(:,~(strcmp(new_names,'POS') | strcmp(new_names,'CHR')));
idx  =intervalsToVector(intervals);
output.N = reshape(sample_sizes(idx), height(output),1);
for i=1:pNum
    output.(['P_PAR', num2str(i)]) = partSignificances(idx, i);
    output.(['CHI_PAR',num2str(i)]) = partScores(idx, i)./(size(PHENO{i},2) .* (1 + partScores(idx, i)./output.N))  ;
end
writetable(output, [save_dir  '/chisq_stats.csv'])
end

function idx = intervalsToVector(intervals)
idx = zeros(intervals(end,2),1);
idx(intervals(:,1)) = 1;
idx = cumsum(idx);
end

%%
function fig = plotPartitionsGWAS(intervals, intStats, chromosome, pThresB, pThres, path)
fig = figure('visible','off');
fig.Position = [100 100 900 900];
hold on
pNum = size(intStats.chisqSignificance, 2);
for i=1:pNum
    sig1 = num2str(sum(intStats.chisqSignificance(:,i)<pThresB));
    sig2 = num2str(sum(intStats.chisqSignificance(:, i)<pThres));
    scatter(intervals(:, 1), -log10(intStats.chisqSignificance(:, i)),'.','DisplayName',['Part. ' num2str(i) ', # significant:', sig1, '(', sig2, ')']);
end
yline(-log10(pThresB), 'DisplayName', 'Bonferroni Threshold');
yline(-log10(pThres), '--', 'DisplayName', '(No correction Threshold)');
title(['Chromosome ' num2str(chromosome) ', Partitions: ' num2str(pNum)])
ylabel('-log10p');
lgd = legend;
set(lgd,'Location','BestOutside');
saveas(fig, [path '_logPlot.svg']);
end


function fig = plotSimpleGWAS(intervals, chisqSignificance, chromosome, pThres, path)
fig = figure;
scatter(intervals(:, 1), -log10(chisqSignificance), 18, '.');
yline(-log10(pThres));

title(['Chromosome ' num2str(chromosome) ', ' num2str(sum(chisqSignificance<pThres)) ' significant intervals out of ' num2str(length(chisqSignificance))])
ylabel('-log10p');
saveas(fig, [path '_logPlot.svg']);
end

function [intSigSnps, sigSnps] = prepareSignificantTablesOnEachPartition(snps, ccaIntStats, ccaIntervals, pThres, save_path)
pNum = size(ccaIntStats.chisqSignificance ,2);
sigSnps= [];
intSigSnps = [];
for i=1:pNum
    partIntStats = struct('chisqSignificance', ccaIntStats.chisqSignificance(:, i));
    [intRet, ret] =  prepareSignificantTables(snps, partIntStats, ccaIntervals, pThres);
    if isempty(ret), continue; end
    intRet.PARTITION = repmat(i,size(intRet, 1),1);
    ret.PARTITION = repmat(i,size(ret, 1),1);
    if size(ret,1) ~= 0
        intSigSnps = [intSigSnps;intRet];
        try
            sigSnps = [sigSnps;ret];
        catch
            disp('h');
        end
    end

end
if ~isempty(sigSnps)
    writetable(sigSnps, [save_path 'significant_snps.csv']);
    writetable(intSigSnps, [save_path 'significant_intervals.csv']);
end
end


function [intSigSnps, sigSnps] = prepareSignificantTables(snps, ccaIntStats, ccaIntervals, pThres, savePath)
if nargin < 5
    savePath = "";
end
significantInts = ccaIntervals(ccaIntStats.chisqSignificance<pThres, :);
starts = significantInts(:,1);
ends = significantInts(:,2);
stSnps = snps(starts,:);
stSnps.Properties.VariableNames = strcat(snps.Properties.VariableNames', '_ST')';
enSnps = snps(ends,:);
enSnps.Properties.VariableNames = strcat(snps.Properties.VariableNames', '_EN')';

% sigSnpsCoef = zeros(size(significantInts, 1),1);
% sigSnps = table();
if isempty(significantInts)
    sigSnps = table();
else
    s = significantInts(:,1);
    e = significantInts(:, 2);
    L = length(s);
    F = cumsum(e-s+1);
    idx = ones(1,F(end));
    idx(1) = s(1);
    idx(1+F(1:L-1)) = s(2:L)-e(1:L-1);
    idx = cumsum(idx);
    sigSnps =  snps(idx, :);
end
intSigSnps = [stSnps, enSnps];
if ~strcmp(savePath, "")
    writetable(sigSnps, [savePath 'significant_snps.csv']);
    writetable(intSigSnps, [savePath 'significant_intervals.csv']);
end
end





