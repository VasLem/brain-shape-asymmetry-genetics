%ENV TO SET:
%%%%%%%%%%%%
%DATA_ROOT: the directory of the DATA (../SAMPLE_DATA/)
%DATASET_INDEX: dataset to use, 1 or 2 (1)
%RESULTS_ROOT: the directory of the results (../results/)
%SCRATCH_ROOT: the directory to use for the intermediate files (../results/)
%THREADS: Number of threads to use (max set by local)
%CHROMOSOME: Chomosome to analyze (21)
%%%%%%%%%%%%
close all;clear;
if ~isdeployed
    restoredefaultpath;
    addpath(genpath('.'));
    addpath(genpath('AIDFUNCTIONS'));
    addpath(genpath('../SNPLIB/'));
    rmpath('../SNPLIB/mexfiles/');% to remove the functions that are causing matlab to crash
    addpath(genpath('SNPLIB-master/mexfiles/'))% where I stored the re-mexed files
end
%%
DATA_DIR = getenv('DATA_ROOT');
if(isempty(DATA_DIR))
    DATA_DIR = '../SAMPLE_DATA/';
end
disp(['Location of data: ', DATA_DIR]);



THREADS = getenv('THREADS');
if(isempty(THREADS))
    THREADS = parcluster('local');
    THREADS = THREADS.NumWorkers;
else
    if ~isnumeric(THREADS)
        THREADS=str2double(THREADS);
    end
end
disp(['Number of threads:', num2str(THREADS)])

CHRS = getenv("CHROMOSOME");
if(isempty(CHRS))
    CHRS = 1:22;
else
    if ~isnumeric(CHRS)
        CHRS=str2double( strsplit(CHRS,','));
    end
end
disp(['Chromosome analyzed:', num2str(CHRS)])

MAX_NUM_FEATS = 0;

DATASET_INDEX = getenv("DATASET_INDEX");
if (isempty(DATASET_INDEX))
    DATASET_INDEX = 1;
else
    disp(DATASET_INDEX)
    if ~isnumeric(DATASET_INDEX)
        DATASET_INDEX = str2double(DATASET_INDEX);
    end
end
disp(['Using dataset:', num2str(DATASET_INDEX)])


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

covGenoPath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];

disp("Loading phenotype and covariates..")
PHENO = load([RESULTS_ROOT, 'hierarchicalClusteringDemo/' DATASET_NAME '/asymmetry_reduction10/ccPriorSegmentation/levels4_mine/phenotype_varThres80.mat']);
load(covGenoPath, 'COV');
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
    PLINK_DATA_OUT = [SCRATCH_CHR_DIR 'plink_data.mat'];
    PLINK_DATA_PROC = ~isdeployed && ~isfile(PLINK_DATA_OUT);
    INT_GENO_OUT = [SCRATCH_CHR_DIR 'intervals_geno.mat'];
    INT_GENO_PROC = ~isdeployed && ~isfile(INT_GENO_OUT);
    META_INT_GENO_OUT = [CHR_DIR 'intervals_info.mat'];
    META_INT_GENO_PROC = ~isfile(META_INT_GENO_OUT);
    CNT_GENO_OUT = [SCRATCH_CHR_DIR 'controlled_geno.mat'];
    CNT_GENO_PROC = ~isdeployed && ~isfile(CNT_GENO_OUT);
    WITH_PART_CCA_OUT = [CHR_DIR 'withPartCCA.mat'];
    WITH_PART_CCA_PROC = ~isfile(WITH_PART_CCA_OUT);
    SAMPLE_SIZES_OUT = [CHR_DIR 'sampleSizes.mat'];
    SAMPLE_SIZES_PROC = ~isfile(SAMPLE_SIZES_OUT);

    if ~isfolder(CHR_DIR), mkdir(CHR_DIR); end
    if ~isfolder(SCRATCH_CHR_DIR), mkdir(SCRATCH_CHR_DIR); end
    
    if ~PLINK_DATA_PROC
        try
            lastwarn('', '');
            load([CHR_DIR 'plink_data_info.mat'], "iid");
            
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
        [snps, samples] = obj.importPLINKDATA([DATA_DIR 'IMAGEN/BRAIN/' UKBIOBANK '/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_' GENO_ID '_chr' num2str(CHR)]);
        geno = obj.UnpackGeno();
        geno(geno==-1) = 255;
        geno = uint8(geno);
        iid = samples.IID;
        toc;
        disp("Saving iids..")
        tic;
        save([CHR_DIR 'plink_data_info.mat'], "iid", '-v7.3')
        clear obj;
        toc;
        if ~isdeployed
            disp("Saving loaded geno..")
            tic;
            save(PLINK_DATA_OUT, 'geno', 'snps', '-v7.3')
            toc;
        end
    end
    %%
    %% Align COV with genotype
    disp("Aligning COV to genotype..");
    tic;
    samplesIId = iid;
    covAssignmentMatrix = (str2double(samplesIId) == str2double(COV.IID)');
    [covGenoIndex, covIndex] = find(covAssignmentMatrix);
    clear covAssignmentMatrix

    covData = COV.DATA(covIndex, :);
    samplesIId = samplesIId(covGenoIndex);

    %%
    genoId = str2double(samplesIId);
    phenoId = str2double(PHENO.preprocPhenoIID);

    toc;
    %% Align phenotype with genotype
    disp("Aligning genotype to phenotype..");
    tic;
    assignmentMatrix = (phenoId == genoId');
    [phenoIndex, genoIndex] = find(assignmentMatrix);
    samplesIId = samplesIId(genoIndex);
    covData = covData(genoIndex,:);
    assert(all(phenoIndex'==1:length(phenoIndex)));
    clear assignmentMatrix
    toc;
    %% Represent SNPs that are not really SNPs (one replaced by many, or many replaced by one)
    % [sortedAf, sortedAfArg] = sort(af);
    % sortedSnps = snps(sortedAfArg, :);
    % fig = figure("Name", "Occurences of InDel including Alleles");
    % hold on
    % yyaxis left
    % lHandle = plot(sortedAf);
    % set(lHandle, {'DisplayName'}, {'Ordered Frequency'});
    %
    % oneToMany = strlength(sortedSnps.(4))>1 & strlength(sortedSnps.(5)) == 1;
    % manyToMany = strlength(sortedSnps.(4))>1 & strlength(sortedSnps.(5)) > 1;
    % manyToOne = strlength(sortedSnps.(4))==1 & strlength(sortedSnps.(5)) > 1;
    % edges = linspace(0, size(manyToOne ,1), 50);
    % oTMH = histcounts(find(oneToMany),edges);
    % mTOH = histcounts(find(manyToOne),edges);
    %
    %
    % sHandle = scatter(find(manyToMany), 0.95 * ones(1, sum(manyToMany))', 'x');
    % sHandle.DisplayName = ['N->N: ', num2str(sum(manyToMany))];
    %
    % yyaxis right
    % bHandle = bar( 0.5*(edges(2:end) + edges(1:end-1)), [oTMH', mTOH'], 1);
    % bHandle(1).FaceColor = 'blue';
    % bHandle(1).DisplayName = ['1->N:  '  num2str(sum(oneToMany))];
    % bHandle(1).EdgeColor = 'k';
    % bHandle(2).FaceColor = 'red';
    % bHandle(2).DisplayName = ['N->1 :', num2str(sum(manyToOne))];
    % bHandle(2).EdgeColor = 'k';
    % xlabel("Sorted Snps")
    % ax = gca;
    % ax.YAxis(2).Color = 'k';
    %
    % legend;
    % hold off
    % savefig(fig, [RESULTS_DIR 'multi_characters_alleles.fig']);
    % saveas(fig, [RESULTS_DIR 'multi_characters_alleles.png']);
    %
    % clear sortedSnps
    %% Remove indels.


    if META_INT_GENO_PROC || INT_GENO_PROC
        if ~exist('geno', 'var')
            disp("Loading geno..")
            tic;
            load([SCRATCH_CHR_DIR 'plink_data.mat'], "geno", "snps");
            toc;
        end
        geno =geno(covGenoIndex, :);
        geno = geno(genoIndex, :);
        disp("Removing SNPs containing InDels..");
        tic;
        indelFilter = strlength(snps.(4))==1 & strlength(snps.(5))==1;
        genoPruned = geno(:, indelFilter);
        snpsPruned = snps(indelFilter, :);
        clear geno
        clear snps

        %% Identify the genetic indices for genes with more than one allele
        [rsids , ~, ic] = unique(snpsPruned.('RSID'));
        alleleCounts = accumarray(ic,1);
        %% Some checks
        %regarding the fact that positions need to be sorted
        assert(all(sort(snpsPruned.POS) == snpsPruned.POS));
        % and that all phenotype ids  correspond to genotype ones
        assert(all((1:length(phenoId)) == phenoIndex'));
        toc;
        %% Partition SNPs
        disp("Computing Intervals..")
        tic;
        intervals = getIntervals(snpsPruned, GENE_SET_METHOD);
        toc;
        disp("Splitting genome into groups, according to intervals information..")
        tic;
        gId = zeros(intervals(end,2),1);
        gId(intervals(:, 1)) = 1;
        gId = cumsum(gId);
        genoInt = splitapply( @(x){x'}, genoPruned', gId );
        toc;
        disp("Saving splitted genome and pruned SNPs...")
        tic;
        save(META_INT_GENO_OUT, 'snpsPruned', 'intervals', 'gId', '-v7.3');
        if ~isdeployed
            save(INT_GENO_OUT, 'genoInt', '-v7.3');
        end
        toc;
        clear genoPruned;
    end

    %%

    if WITH_PART_CCA_PROC
        if CNT_GENO_PROC
            if ~exist('genoInt', 'var')
                disp("Loading genoInt..")
                tic;
                load(INT_GENO_OUT, 'genoInt');
                toc;
            end
            disp("Controlling genome for COV based on intervals..")
            tic;
            genoControlledInt = controlGenoCovariates(genoInt, covData);
            toc;
            if ~isdeployed
                disp("Saving controlled genome..")
                tic;
                save(CNT_GENO_OUT, 'genoControlledInt', '-v7.3');
                toc;
            end
            if ~SAMPLE_SIZES_PROC
                clear genoInt
            end
        end
    
    end
    if SAMPLE_SIZES_PROC
        if ~exist('genoInt', 'var')
                disp("Loading genoInt..")
                tic;
                load(INT_GENO_OUT, 'genoInt');
                toc;
        end
        sampleSizes = getSampleSizes(genoInt);
        save(SAMPLE_SIZES_OUT,'sampleSizes', '-v7.3')
    end
    %%
    if ~exist('intervals', 'var')
        disp("Loading intervals from memory..")
        tic;
        load(META_INT_GENO_OUT,  'intervals', 'gId');
        toc;
    end
    gTLPartsPThres = NO_PARTITION_THRES / length(PHENO.clusterPCAPhenoFeatures);
    if ~WITH_PART_CCA_PROC
        load(WITH_PART_CCA_OUT,  'gTLPartStats', 'gTLPartIntStats');
        disp("Loaded Computed CCA with phenotypic partinioning");
    else
        if ~exist('genoControlledInt', 'var')
            disp("Loading genoControlledInt..")
            tic;
            load(CNT_GENO_OUT, 'genoControlledInt');
            toc;
        end
        disp("Computing CCA with phenotypic partitioning..")
        tic;
        [gTLPartStats, gTLPartIntStats] = runCCAOnEachPartition(PHENO, genoControlledInt, intervals, gId, SCRATCH_CHR_DIR, MAX_NUM_FEATS);
        toc;
        disp("Saving CCA results..")
        tic;
        save(WITH_PART_CCA_OUT, 'gTLPartStats', 'gTLPartIntStats', '-v7.3');
        WITH_PART_CCA_PROC = false;
        toc;
    end
    %%
    disp("Plotting results..")
    tic;
    plotSimpleGWAS(intervals, gTLPartIntStats.chiSqSignificance(1, :), CHR, NO_PARTITION_THRES,  [CHR_DIR  'noPartition_feats' num2str(MAX_NUM_FEATS)]);
    plotPartitionsGWAS(intervals, gTLPartIntStats, CHR, gTLPartsPThres, [CHR_DIR 'PartitionedGTL_feats' num2str(MAX_NUM_FEATS)]);
    plotPartitionsGWAS(intervals, gTLPartIntStats, CHR, NO_PARTITION_THRES, [CHR_DIR 'PartitionedGTL_NOBF_feats' num2str(MAX_NUM_FEATS)]);
    toc;
    %% Significant SNPs tables extraction
    
    %%
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
        disp("Loading snpsPruned..")
        tic;
        load(SAMPLE_SIZES_OUT);
        toc;
    end
    tic;
    saveLDRegressionTablesOnEachPartition(snpsPruned,sampleSizes, gTLPartIntStats.chiSqSignificance, intervals, CHR_DIR);
    tic;
    %
    disp("End of computation.")
    %%
    if ~isdeployed
        close all
    end
    toc;
end
function output = saveLDRegressionTablesOnEachPartition(snpsPruned, sample_sizes, partSignificances, intervals, save_dir)
pNum = size(partSignificances ,1);
output = removevars(renamevars(snpsPruned,["ALT", "REF"],["A2", "A1"]), ["POS", "CHR"]);
idx  =intervalsToVector(intervals);
output.N = sample_sizes(idx);
output.SIGN = ones(height(snpsPruned),1);
for i=1:pNum
    output.(['P_PAR', num2str(i)]) = partSignificances(i, idx)';
end
writetable(output, [save_dir  '/chisq_stats.csv'])
end

function idx = intervalsToVector(intervals)
    idx = zeros(intervals(end,2),1);
    idx(intervals(:,1)) = 1;
    idx = cumsum(idx);
end

%%
function fig = plotPartitionsGWAS(intervals, intStats, chromosome, pThres, path)
fig = figure;
fig.Position = [100 100 900 900];
hold on
pNum = size(intStats.chiSqSignificance, 1);
for i=1:pNum
    scatter(intervals(:, 1), -log10(intStats.chiSqSignificance(i, :)),'.','DisplayName',['Part. ' num2str(i) ', # significant:', num2str(sum(intStats.chiSqSignificance(i, :)<pThres))]);
end
yline(-log10(pThres), 'DisplayName', 'Threshold');

title(['Chromosome ' num2str(chromosome) ', Partitions: ' num2str(pNum)])
ylabel('-log10p');

lgd = legend;
set(lgd,'Location','BestOutside');

%%
saveas(fig, [path '_logPlot.png']);
end


function fig = plotSimpleGWAS(intervals, chiSqSignificance, chromosome, pThres, path)
%%
fig = figure;
scatter(intervals(:, 1), -log10(chiSqSignificance),'.');
yline(-log10(pThres));

title(['Chromosome ' num2str(chromosome) ', ' num2str(sum(chiSqSignificance<pThres)) ' significant intervals out of ' num2str(length(chiSqSignificance))])
ylabel('-log10p');

%%
saveas(fig, [path '_logPlot.png']);
end

function [intSigSnps, sigSnps] = prepareSignificantTablesOnEachPartition(snps, ccaIntStats, ccaIntervals, pThres, save_path)
pNum = size(ccaIntStats.chiSqSignificance ,1);
sigSnps= [];
intSigSnps = [];
for i=1:pNum
    partIntStats = struct('chiSqSignificance', ccaIntStats.chiSqSignificance(i, :));
    %     partStats.coeffs = ccaStats.coeffs(i, :);
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

function [stats, intStats, intervals] = runCCAOnEachPartition(PHENO, geno, intervals, intIdVec, scratch_dir, maxNumPhenoFeats)
if nargin < 6
    maxNumPhenoFeats = 0;
end
pnum = length(PHENO.clusterPCAPhenoFeatures);
gnum = intervals(end,2);
inum = size(intervals, 1);
chisq = zeros(pnum, gnum);
chisqSig = zeros(pnum, gnum);
intChisqSig = zeros(pnum, inum);
intChisq = zeros(pnum, inum);
for k=pnum:-1:1
    phenoPart = PHENO.clusterPCAPhenoFeatures{k};
    if maxNumPhenoFeats ~= 0
        phenoPart = phenoPart(:,1:min(size(phenoPart,2),maxNumPhenoFeats));
    end
    disp(['Partition ', num2str(k), ', number of components: ', num2str(size(phenoPart,2))]);
    PART_DIR = [scratch_dir 'par' num2str(k) '/'];
    if ~isfolder(PART_DIR), mkdir(PART_DIR); end
    try
        load([PART_DIR 'data.mat'], 's', 'i');
        disp("Loaded from disk")
    catch
        tic;
        [s, i] = runCCA(phenoPart, geno, intervals, intIdVec);
        toc;

        save([PART_DIR 'data.mat'], 's', 'i', '-v7.3');
    end

    chisqSig(k, :) = s.chiSqSignificance;
    chisq(k, :) = s.chisq;
    intChisqSig(k, :) = i.chiSqSignificance;
    intChisq(k, :) = i.chisq;
end

intStats.chiSqSignificance = intChisSq;
intStats.chisq = intChisq;
stats.chiSqSignificance = chisq;
stats.chisq = chisq;
end





function controlledGenoPart = controlGenoCovariates(genoInt, COV)
controlledGenoPart = genoInt;
ppb = ParforProgressbar(length(genoInt));
parfor i=1:length(genoInt)
    controlledGenoPart{i} = single(nan * zeros(size(genoInt{i})));
    selection = ~any(genoInt{i} == 255,2);
    controlledGenoPart{i}(selection, :) = single(getResiduals(COV(selection, :),  single(genoInt{i}(selection, :))));
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
    sampleSizes = cellfun(@(x)(min(sum(isfinite(x)))), geno);
end
