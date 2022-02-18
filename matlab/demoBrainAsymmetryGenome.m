%ENV TO SET:
%%%%%%%%%%%%
%DATA_ROOT: the directory of the DATA (../SAMPLE_DATA/)
%DATASET_INDEX: dataset to use, 1 or 2 (1)
%RESULTS_ROOT: the directory of the results (../results/)
%SCRATCH_ROOT: the directory to use for the intermediate files (../results/)
%THREADS: Number of threads to use (max set by local)
%CHROMOSOME: Chomosome to analyze (21)
%BLOCK_SIZE: Block Size for CCA (2000 when MEDIAN_IMPUTE is on else 20000)
%MEDIAN_IMPUTE: Whether to perform median imputation (very fast) or not (very slow) (0)
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
MAX_NUM_FEATS = 0;
NO_PARTITION_THRES = 5*10^-8; % European in LD score
DEFAULT_CHRS = 1:22;
DEFAULT_DATASET_INDEX = 1;
DEFAULT_MEDIAN_IMPUTE = 1;

MEDIAN_IMPUTE = getenv('MEDIAN_IMPUTE');
if(isempty(MEDIAN_IMPUTE))
    MEDIAN_IMPUTE = DEFAULT_MEDIAN_IMPUTE;
end

if MEDIAN_IMPUTE == 0
    DEFAULT_BLOCK_SIZE= 20000;
else
    DEFAULT_BLOCK_SIZE = 2000;
end

BLOCK_SIZE = getenv('BLOCK_SIZE');
if(isempty(BLOCK_SIZE))
    BLOCK_SIZE = DEFAULT_BLOCK_SIZE;
end


DATA_DIR = getenv('DATA_ROOT');
if(isempty(DATA_DIR))
    DATA_DIR = '../SAMPLE_DATA/';
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

disp(['Number of threads:', num2str(THREADS)])
disp(['Location of data: ', DATA_DIR])
disp(['Median Imputation: ', num2str(MEDIAN_IMPUTE)])
disp(['Block Size: ', num2str(BLOCK_SIZE)]);
disp(['Using dataset:', num2str(DATASET_INDEX)])
disp(['Chromosome analyzed:', num2str(CHRS)])

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
RESULTS_DIR = [RESULTS_ROOT, 'genomeDemo/' DATASET_NAME '/'];
disp(['Location of results: ', RESULTS_DIR]);
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end

SCRATCH_ROOT = getenv('SCRATCH_ROOT');
if(isempty(SCRATCH_ROOT))
    SCRATCH_ROOT = '../results/';
end

disp(['Location of intermediate files: ', SCRATCH_ROOT]);
SCRATCH_DIR = [SCRATCH_ROOT, '/genomeDemo/' DATASET_NAME '/'];
if ~isfolder(SCRATCH_DIR), mkdir(SCRATCH_DIR); end

COV_GENO_PATH = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];
disp("Loading phenotype..")
PHENO_PATH = [RESULTS_ROOT, 'hierarchicalClusteringDemo/' DATASET_NAME '/asymmetry_reduction10/ccPriorSegmentation/levels4_mine/phenotype_varThres80.mat'];
PHENO = load(PHENO_PATH);
%%
disp("Loading covariates..")
load(COV_GENO_PATH, 'COV');
disp("Initializing genomic analysis..")
try
    parpool(THREADS);
catch
end
for CHR_IND=1:length(CHRS)
    tic;
    clearvars -except -regexp ^[A-Z|_]+$
    CHR = CHRS(CHR_IND);
    disp(['CHR:' , num2str(CHR)]);
    GENE_SET_METHOD = 'perSNP';
    CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
    SCRATCH_CHR_DIR = [SCRATCH_DIR 'chr' num2str(CHR) '/'];

    PLINK_DATA_INFO_OUT = [CHR_DIR 'plink_data_info.mat'];
    PLINK_DATA_PROC = isdeployed || ~isfile(PLINK_DATA_INFO_OUT);
    META_INT_GENO_OUT = [CHR_DIR 'intervals_info.mat'];
    META_INT_GENO_PROC = ~isfile(META_INT_GENO_OUT);
    WITH_PART_CCA_OUT = [CHR_DIR 'withPartCCA.mat'];
    WITH_PART_CCA_PROC = ~isfile(WITH_PART_CCA_OUT);
    SAMPLE_SIZES_OUT = [CHR_DIR 'sampleSizes.mat'];

    if ~isfolder(CHR_DIR), mkdir(CHR_DIR); end
    if ~isfolder(SCRATCH_CHR_DIR), mkdir(SCRATCH_CHR_DIR); end
    PLINK_DATA_PROC = PLINK_DATA_PROC || META_INT_GENO_PROC || WITH_PART_CCA_PROC;
    META_INT_GENO_PROC =  META_INT_GENO_PROC || ~isfile(PLINK_DATA_INFO_OUT);
    WITH_PART_CCA_PROC = WITH_PART_CCA_PROC || META_INT_GENO_PROC || ~isfile(PLINK_DATA_INFO_OUT);


    if ~PLINK_DATA_PROC
        try
            lastwarn('', '');
            load([CHR_DIR 'plink_data_info.mat'], "iid", 'phenoIndex');

            [warnMsg, warnId] = lastwarn();
            if ~isempty(warnId)
                error(warnMsg, warnId);
            end
            disp("Genomic Data loaded from disk")
        catch
            PLINK_DATA_PROC = true;
        end
    end
    if PLINK_DATA_PROC
        disp("Loading PLINK preprocessed data..")
        tic;
        obj = SNPLIB();
        if ~isdeployed
            obj.nThreads = THREADS;
        end
        GENO_PATH = [DATA_DIR 'IMAGEN/BRAIN/' UKBIOBANK '/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_' GENO_ID '_chr' num2str(CHR)];
        [snps, samples] = obj.importPLINKDATA(GENO_PATH);
        geno = obj.UnpackGeno();
       

        clear obj;

        f = figure();
        histogram(100 * (sum(geno==255) / size(geno, 1) ));
        xlabel('Missing Values (%)')
        ylabel("SNPs #")
        saveas(f,[CHR_DIR,'missing_histogram.svg'])
        toc;
        % Align phenotype with genotype
        disp("Aligning genotype to phenotype..");
        tic;
        iid = samples.IID;
        genoId = str2double(iid);
        phenoId = str2double(PHENO.preprocPhenoIID);
        assignmentMatrix = (phenoId == genoId');
        [phenoIndex, genoIndex] = find(assignmentMatrix);
        iid = iid(genoIndex);
        geno = geno(genoIndex, :);
        assert(all(phenoIndex'==1:length(phenoIndex)));
        clear assignmentMatrix
        toc;

        % Remove indels.
        disp("Removing SNPs containing InDels..");
        tic;
        indelFilter = strlength(snps.(4))==1 & strlength(snps.(5))==1;
        genoPruned = geno(:, indelFilter);
        snpsPruned = snps(indelFilter, :);
        clear geno
        clear snps

        % Some checks
        %regarding the fact that positions need to be sorted
        assert(all(sort(snpsPruned.POS) == snpsPruned.POS));
        % and that all phenotype ids  correspond to genotype ones
        assert(all((1:length(phenoId)) == phenoIndex'));
        toc;

        % Partition SNPs
        disp("Computing Intervals..")
        tic;
        intervals = getIntervals(snpsPruned, GENE_SET_METHOD);
        toc;
        %Align COV with genotype
        disp("Aligning COV to genotype and phenotype..");
        tic;
        samplesIId = iid;
        covAssignmentMatrix = (str2double(samplesIId) == str2double(COV.IID)');
        [covGenoIndex, covIndex] = find(covAssignmentMatrix);
        clear covAssignmentMatrix
        iid = samplesIId(covGenoIndex);
        phenoIndex = phenoIndex(covGenoIndex);

        genoPruned = genoPruned(covGenoIndex, :);
        covData = COV.DATA(covIndex);
        toc;
        if ~isdeployed
            disp("Saving aligned iids..")
            tic;
            save(PLINK_DATA_INFO_OUT, "iid", "phenoIndex", '-v7.3')
        end
        toc;
        if MEDIAN_IMPUTE
            disp("Imputing missing using median..")
            tic;
            parfor i=1:size(genoPruned,2)
                s = genoPruned(:, i);
                m =s~=255;
                medians(i) = round(median(s(m)));
                s(~m) = medians(i);
                genoPruned(:, i) = s;
            end
            toc;
        end

        disp("Splitting genome into groups, according to intervals information..")
        tic;
        gId = zeros(intervals(end,2),1);
        gId(intervals(:, 1)) = 1;
        gId = cumsum(gId);
        genoPruned = splitapply( @(x){x'}, genoPruned', gId );

        toc;
        if ~isdeployed
            disp("Saving pruned SNPs...")
            tic;
            save(META_INT_GENO_OUT, 'snpsPruned', 'intervals', 'gId', '-v7.3');
            toc;
        end
        disp("Computing sample sizes per SNP..")
        sampleSizes = getSampleSizes(genoPruned);
        if ~isdeployed
            save(SAMPLE_SIZES_OUT,'sampleSizes', '-v7.3')
        end

        %         disp("Controlling genome for covariates based on intervals..")
        %         tic;
        %         genoInt = controlGenoCovariates(genoInt, covData);
        %         toc;
    end
    %%
    if ~exist('intervals', 'var')
        disp("Loading intervals from memory..")
        tic;
        load(META_INT_GENO_OUT,  'intervals', 'gId');
        toc;
    end
    pheno = getPheno(PHENO, phenoIndex);
    gTLPartsPThres = NO_PARTITION_THRES / length(pheno);
    %%
    if ~WITH_PART_CCA_PROC
        load(WITH_PART_CCA_OUT,  'gTLPartStats', 'gTLPartIntStats');
        disp("Loaded Computed CCA with phenotypic partinioning");
    else
         %%
        disp("Computing CCA with phenotypic partitioning..")
        tic;
       
        [gTLPartStats, gTLPartIntStats] = runCCA(pheno,  genoPruned, intervals, gId, BLOCK_SIZE);
        toc;
        disp("Saving CCA results..")
        tic;
        save(WITH_PART_CCA_OUT, 'gTLPartStats', 'gTLPartIntStats', '-v7.3');
        WITH_PART_CCA_PROC = false;
        toc;
        clear genoPruned
    end

    %%
    if UPDATE_FIGS
        disp("Plotting results..")
        tic;
        plotSimpleGWAS(intervals, gTLPartIntStats.chiSqSignificance(:, 1), CHR, NO_PARTITION_THRES,  [CHR_DIR  'noPartition_feats' num2str(MAX_NUM_FEATS)]);
        plotPartitionsGWAS(intervals, gTLPartIntStats, CHR, gTLPartsPThres, NO_PARTITION_THRES, [CHR_DIR 'PartitionedGTL_feats' num2str(MAX_NUM_FEATS)]);
        toc;
    end
    %% Significant SNPs tables extraction
    if ~exist('snpsPruned', 'var')
        disp("Loading snpsPruned..")
        tic;
        load(META_INT_GENO_OUT, 'snpsPruned');
        toc;
    end
    disp("Extracting significant SNPs tables..")
    tic;
    prepareSignificantTablesOnEachPartition(snpsPruned,  gTLPartIntStats , intervals, gTLPartsPThres, [CHR_DIR, 'PartitionedGTLWithBC_feats' num2str(MAX_NUM_FEATS)]);
    prepareSignificantTablesOnEachPartition(snpsPruned,  gTLPartIntStats , intervals, NO_PARTITION_THRES, [CHR_DIR, 'PartitionedGTLWoutBC_feats' num2str(MAX_NUM_FEATS)]);
    toc;
    %% LD regression csv
    disp("Extracting partition specific signigicant SNPs tables..")
    if ~exist('sampleSizes', 'var')
        disp("Loading sampleSizes..")
        tic;
        load(SAMPLE_SIZES_OUT);
        toc;
    end
    tic;
    saveLDRegressionTablesOnEachPartition(snpsPruned, pheno, sampleSizes, gTLPartIntStats.chisq, gTLPartIntStats.chiSqSignificance, intervals, CHR_DIR);
    tic;
    %
    disp("End of computation.")
    %%
    if ~isdeployed
        close all
    end
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
output.N = sample_sizes(idx);
output.SIGN = ones(height(snpsPruned),1);
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
pNum = size(intStats.chiSqSignificance, 2);
for i=1:pNum
    sig1 = num2str(sum(intStats.chiSqSignificance(:,i)<pThresB));
    sig2 = num2str(sum(intStats.chiSqSignificance(:, i)<pThres));
    scatter(intervals(:, 1), -log10(intStats.chiSqSignificance(:, i)),'.','DisplayName',['Part. ' num2str(i) ', # significant:', sig1, '(', sig2, ')']);
end
yline(-log10(pThresB), 'DisplayName', 'Bonferroni Threshold');
yline(-log10(pThres), '--', 'DisplayName', '(No correction Threshold)');
title(['Chromosome ' num2str(chromosome) ', Partitions: ' num2str(pNum)])
ylabel('-log10p');
lgd = legend;
set(lgd,'Location','BestOutside');
saveas(fig, [path '_logPlot.svg']);
end


function fig = plotSimpleGWAS(intervals, chiSqSignificance, chromosome, pThres, path)
fig = figure;
scatter(intervals(:, 1), -log10(chiSqSignificance), 18, '.');
yline(-log10(pThres));

title(['Chromosome ' num2str(chromosome) ', ' num2str(sum(chiSqSignificance<pThres)) ' significant intervals out of ' num2str(length(chiSqSignificance))])
ylabel('-log10p');
saveas(fig, [path '_logPlot.svg']);
end

function [intSigSnps, sigSnps] = prepareSignificantTablesOnEachPartition(snps, ccaIntStats, ccaIntervals, pThres, save_path)
pNum = size(ccaIntStats.chiSqSignificance ,2);
sigSnps= [];
intSigSnps = [];
for i=1:pNum
    partIntStats = struct('chiSqSignificance', ccaIntStats.chiSqSignificance(:, i));
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
significantInts = ccaIntervals(ccaIntStats.chiSqSignificance<pThres, :);
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
% for c=1:size(significantInts,1)
%     intSnps = snps(starts(c):ends(c),:);
%     intCoeffs = ccaStats.coeffs(starts(c):ends(c));
%     [~, maxCoefArg] = max(abs(intCoeffs));
%     sigSnpsCoef(c) = intCoeffs(maxCoefArg);
%     sigSnps(c, :) = intSnps(maxCoefArg, :);
% end
% sigSnps.COEF = sigSnpsCoef;
intSigSnps = [stSnps, enSnps];
if ~strcmp(savePath, "")
    writetable(sigSnps, [savePath 'significant_snps.csv']);
    writetable(intSigSnps, [savePath 'significant_intervals.csv']);
end
end

function ret = getPheno(PHENO, phenoIndex)
    pnum = length(PHENO.clusterPCAPhenoFeatures);
    ret = cell(length(pnum),1);
    for k=1:pnum
        phenoPart = PHENO.clusterPCAPhenoFeatures{k};
        phenoPart = phenoPart(phenoIndex, :);
        ret{k} = phenoPart;
    end
end


function genoInt = controlGenoCovariates(genoInt, COV)
ppb = ParforProgressbar(length(genoInt));
parfor i=1:length(genoInt)
    selection = ~any(genoInt{i} == 255,2);
    genoInt{i}(selection, :) = single(getResiduals(COV(selection, :),  single(genoInt{i}(selection, :))));
    ppb.increment();
end
delete(ppb);
end

function intervals = getIntervals(snps, geneSetMethod, namedArgs)
% geneSetMethod:
%  slidingWindow: overlapping set of SNPS are taken based on the provided windowSize
%  perSnp: multiallelic SNPS are considered in the same CCA
if nargin < 2
    geneSetMethod = 'perSNP';
end
if nargin < 3
    namedArgs.windowSize = 20000;
end
switch geneSetMethod
    case 'slidingWindow'
        snpCnt = 1;
        cnt = 1;
        intSlidWindow = zeros(size(snps, 1),2);
        while snpCnt<size(snps, 1)
            startInd = snpCnt;
            startPos = snps.POS(snpCnt);
            offset = find(snps.POS(snpCnt:end) > startPos + namedArgs.windowSize, 1,  'first');
            if size(offset)==0, offset = size(snps, 1) + 1; end
            endInd = snpCnt + offset - 1;
            snpCnt = endInd + 1;
            intSlidWindow(cnt, :) = [startInd, endInd];
            cnt = cnt +1;
        end
        intervals = intSlidWindow(1:cnt-1, 2);
    case 'perSNP'
        [uniquePos, ia, ~] = unique(snps.POS);
        intervals = zeros(length(uniquePos), 2);
        intervals(:, 1) = ia(1:end);
        intervals(:, 2) = [ia(2:end)-1;length(snps.POS)];
end
end

function sampleSizes = getSampleSizes(geno)
sampleSizes = cellfun(@(x)(min(sum(x~=255))), geno);
end
