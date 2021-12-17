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

CHR = getenv("CHROMOSOME");
if(isempty(CHR))
    CHR = 21;
else
    if ~isnumeric(CHR)
        CHR=str2double(CHR);
    end
end
disp(['Chromosome analyzed:', num2str(CHR)])

MAX_NUM_FEATS = 0;

DATASET_INDEX = getenv("DATASET_INDEX");
if (isempty(DATASET_INDEX))
    DATASET_INDEX = 2;
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
pheno = load([RESULTS_ROOT, 'hierarchicalClusteringDemo/' DATASET_NAME '/asymmetry_reduction10/ccPriorSegmentation/levels4_mine/phenotype_varThres80.mat']);
load(covGenoPath, 'COV');
covariates = COV;
disp("Initializing genomic analysis..")
try
    parpool(THREADS);
catch
end

disp(['CHR:' , num2str(CHR)]);
GENE_SET_METHOD = 'perSNP';
CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
SCRATCH_CHR_DIR = [SCRATCH_DIR 'chr' num2str(CHR) '/'];
INT_GENO_OUT = [SCRATCH_CHR_DIR 'intervals_geno.mat'];
INT_GENO_PROC = ~isfile(INT_GENO_OUT);
META_INT_GENO_OUT = [CHR_DIR 'intervals_info.mat'];
META_INT_GENO_PROC = ~isfile(META_INT_GENO_OUT);
CNT_GENO_OUT = [SCRATCH_CHR_DIR 'controlled_geno.mat'];
CNT_GENO_PROC = ~isfile(CNT_GENO_OUT);
NO_PART_CCA_OUT = [CHR_DIR 'noPartCCA.mat'];
NO_PART_CCA_PROC = ~isfile(NO_PART_CCA_OUT);
WITH_PART_CCA_OUT = [CHR_DIR 'withPartCCA.mat'];
WITH_PART_CCA_PROC = ~isfile(WITH_PART_CCA_OUT);

if ~isfolder(CHR_DIR), mkdir(CHR_DIR); end
if ~isfolder(SCRATCH_CHR_DIR), mkdir(SCRATCH_CHR_DIR); end
try
    lastwarn('', '');
    load([CHR_DIR 'plink_data_info.mat'], "iid");
    samples.IID = iid;
    if META_INT_GENO_PROC
        load([SCRATCH_CHR_DIR 'plink_data.mat'], "geno", "snps", "samples");
    end
    [warnMsg, warnId] = lastwarn();
    if ~isempty(warnId)
        error(warnMsg, warnId);
    end
    disp("Genomic Data loaded from disk")
catch
    disp("Loading PLINK preprocessed data..")
    obj = SNPLIB();
    if ~isdeployed
        obj.nThreads = THREADS;
    end
    [snps, samples] = obj.importPLINKDATA([DATA_DIR 'IMAGEN/BRAIN/' UKBIOBANK '/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_' GENO_ID '_chr' num2str(CHR)]);
    geno = obj.UnpackGeno();
    geno(geno==-1) = 255;
    geno = uint8(geno);
    iid = samples.IID;
    if ~isdeployed
        disp("Computing LD scores..")
        ld = obj.CalcLDScores(snps);
        disp("Computing Allele Frequencies..")
        af = obj.CalcAlleleFrequency();
        save([CHR_DIR 'plink_data_info.mat'], 'af', 'ld',"iid", '-v7.3')
    else
     save([CHR_DIR 'plink_data_info.mat'], "iid", '-v7.3')
    end
    disp("Saving to disk..")
    clear obj;
    save([SCRATCH_CHR_DIR 'plink_data.mat'], 'geno', 'snps', 'samples', '-v7.3')
    INT_GENO_PROC = 1;
    CNT_GENO_PROC = 1;
    NO_PART_CCA_PROC = 1;
    WITH_PART_CCA_PROC = 1;
end
%%
%% Align covariates with genotype
disp("Aligning covariates to genotype..");
covAssignmentMatrix = (str2double(samples.IID) == str2double(covariates.IID)');
[covGenoIndex, covIndex] = find(covAssignmentMatrix);
clear covAssignmentMatrix
covData = covariates.DATA(covIndex, :);
if META_INT_GENO_PROC
    geno =geno(covGenoIndex, :);
    samples = samples(covGenoIndex, :);
end


%%
genoId = str2double(samples.IID);
phenoId = str2double(pheno.preprocPhenoIID);


%% Align phenotype with genotype
disp("Aligning genotype to phenotype..");
assignmentMatrix = (phenoId == genoId');
[phenoIndex, genoIndex] = find(assignmentMatrix);
if META_INT_GENO_PROC
    geno = geno(genoIndex, :);
    samples = samples(genoIndex, :);
end
covData = covData(genoIndex,:);
assert(all(phenoIndex'==1:length(phenoIndex)));
clear assignmentMatrix

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
disp("Removing SNPs containing InDels..");
indelFilter = strlength(snps.(4))==1 & strlength(snps.(5))==1;
if META_INT_GENO_PROC
    genoPruned = geno(:, indelFilter);
    clear geno
    snpsPruned = snps(indelFilter, :);
    clear snps
end


%% Identify the genetic indices for genes with more than one allele
[rsids , ~, ic] = unique(snpsPruned.('RSID'));
alleleCounts = accumarray(ic,1);
%% Some checks
%regarding the fact that positions need to be sorted
assert(all(sort(snpsPruned.POS) == snpsPruned.POS));
% and that all phenotype ids  correspond to genotype ones
assert(all((1:length(phenoId)) == phenoIndex'));
%% Partition SNPs


if ~META_INT_GENO_PROC
    if NO_PART_CCA_PROC || WITH_PART_CCA_PROC
        if CNT_GENO_PROC
            load(META_INT_GENO_OUT, 'snpsPruned', 'intervals', 'gId');
            load(INT_GENO_OUT, 'genoInt');
        else
            load(META_INT_GENO_OUT, 'snpsPruned', 'intervals', 'gId');
        end

    else
        load(META_INT_GENO_OUT, 'intervals', 'gId');
    end
    disp("Loaded intervals.");
else
    disp("Computing Intervals..")
    intervals = getIntervals(snpsPruned, GENE_SET_METHOD);
    disp("Splitting genome into groups, according to intervals information..")
    gId = zeros(intervals(end,2),1);
    gId(intervals(:, 1)) = 1;
    gId = cumsum(gId);
    genoInt = splitapply( @(x){x'}, genoPruned', gId );
    save(META_INT_GENO_OUT, 'snpsPruned', 'intervals', 'gId', '-v7.3');
    save(INT_GENO_OUT, 'genoInt', '-v7.3');
    clear genoPruned;
end

%%

if NO_PART_CCA_PROC || WITH_PART_CCA_PROC
    if ~CNT_GENO_PROC
        load(CNT_GENO_OUT, 'genoControlledInt');
        disp("Loaded controlled genome for covariates based on intervals.")
    else
        disp("Controlling genome for covariates based on intervals..")
        genoControlledInt = controlGenoCovariates(genoInt, covData);
        save(CNT_GENO_OUT, 'genoControlledInt', '-v7.3');
        clear genoInt
    end

end

%%

if ~NO_PART_CCA_PROC
    load(NO_PART_CCA_OUT, 'noPartitionStats', 'noPartitionIntStats');
    disp("Loaded Computed CCA without phenotypic partitioning");
else
    disp("Computing CCA without phenotypic partitioning..")
    % Without Taking Partitioning Into Consideration
    noPartPheno = pheno.clusterPCAPhenoFeatures{1};
    if MAX_NUM_FEATS~=0
        noPartPheno = noPartPheno(:,1:MAX_NUM_FEATS);
    end
    [noPartitionStats, noPartitionIntStats] = runCCA(noPartPheno, genoControlledInt, intervals, gId);
    save(NO_PART_CCA_OUT, 'noPartitionStats', 'noPartitionIntStats', '-v7.3');
    clear noPartPheno;
    NO_PART_CCA_PROC = False;
end

plotSimpleGWAS(intervals, noPartitionIntStats, CHR, NO_PARTITION_THRES,  [CHR_DIR  'noPartition_feats' num2str(MAX_NUM_FEATS)]);
%%
gTLPartsPThres = NO_PARTITION_THRES / length(pheno.clusterPCAPhenoFeatures);
if ~WITH_PART_CCA_PROC
    load(WITH_PART_CCA_OUT,  'gTLPartStats', 'gTLPartIntStats');
    disp("Loaded Computed CCA with phenotypic partinioning");
else

    % With Global-To-Local Partitioning Into Consideration
    disp("Computing CCA with phenotypic partitioning..")
    [gTLPartStats, gTLPartIntStats] = runCCAOnEachPartition(pheno, genoControlledInt, intervals, gId, SCRATCH_CHR_DIR, MAX_NUM_FEATS);


    save(WITH_PART_CCA_OUT, 'gTLPartStats', 'gTLPartIntStats', '-v7.3');
    WITH_PART_CCA_PROC = False;
end
if ~NO_PART_CCA_PROC && ~WITH_PART_CCA_PROC
    delete(CNT_GENO_OUT);
end
%%
plotPartitionsGWAS(intervals, gTLPartIntStats, CHR, gTLPartsPThres, [CHR_DIR 'PartitionedGTL_feats' num2str(MAX_NUM_FEATS)]);
plotPartitionsGWAS(intervals, gTLPartIntStats, CHR, NO_PARTITION_THRES, [CHR_DIR 'PartitionedGTL_NOBF_feats' num2str(MAX_NUM_FEATS)]);

%%
%% Significant SNPs tables extraction
disp("Extracting significant SNPs tables..")
%%
prepareSignificantTablesOnEachPartition(snpsPruned,  gTLPartIntStats , intervals, gTLPartsPThres, [CHR_DIR, 'PartitionedGTLWithBC_feats' num2str(MAX_NUM_FEATS)]);
%%
prepareSignificantTablesOnEachPartition(snpsPruned,  gTLPartIntStats , intervals, NO_PARTITION_THRES, [CHR_DIR, 'PartitionedGTLWoutBC_feats' num2str(MAX_NUM_FEATS)]);
%
disp("End of computation.")
%%
if ~isdeployed
    close all
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


function fig = plotSimpleGWAS(intervals, intStats, chromosome, pThres, path)
%%
fig = figure;
scatter(intervals(:, 1), -log10(intStats.chiSqSignificance),'.');
yline(-log10(pThres));

title(['Chromosome ' num2str(chromosome) ', ' num2str(sum(intStats.chiSqSignificance<pThres)) ' significant intervals out of ' num2str(length(intStats.chiSqSignificance))])
ylabel('-log10p');

%%
saveas(fig, [path '_logPlot.png']);
end

function [intSigSnps, sigSnps] = prepareSignificantTablesOnEachPartition(snps, ccaIntStats, ccaIntervals, pThres, save_path)
pNum = size(ccaIntStats.chiSqSignificance ,1);
sigSnps= [];
intSigSnps = [];
parfor i=1:pNum
    partIntStats = struct('chiSqSignificance',ccaIntStats.chiSqSignificance(i, :));
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
sigSnps = table();
for c=1:size(significantInts,1)
    intSnps = snps(starts(c):ends(c),:);
    %     intCoeffs = ccaStats.coeffs(starts(c):ends(c));
    %     [~, maxCoefArg] = max(abs(intCoeffs));
    %     sigSnpsCoef(c) = intCoeffs(maxCoefArg);
    %     sigSnps(c, :) = intSnps(maxCoefArg, :);
    sigSnps(c, :) = intSnps(1, :);
end
% sigSnps.COEF = sigSnpsCoef;
intSigSnps = [stSnps, enSnps];
if ~strcmp(savePath, "")
    writetable(sigSnps, [savePath 'significant_snps.csv']);
    writetable(intSigSnps, [savePath 'significant_intervals.csv']);
end
end

function [stats, intStats, intervals] = runCCAOnEachPartition(pheno, geno, intervals, intIdVec, scratch_dir, maxNumPhenoFeats)
if nargin < 6
    maxNumPhenoFeats = 0;
end
pnum = length(pheno.clusterPCAPhenoFeatures);
gnum = intervals(end,2);
inum = size(intervals, 1);
chisq = zeros(pnum, gnum);
% coeffs = zeros(pnum, gnum);
intChisSq = zeros(pnum, inum);
for k=pnum:-1:1
    phenoPart = pheno.clusterPCAPhenoFeatures{k};
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

    chisq(k, :) = s.chiSqSignificance;
    %     coeffs(k, :) = s.coeffs;
    intChisSq(k, :) = i.chiSqSignificance;
    %     intCoeffs{k} = i.coeffs;
end

intStats.chiSqSignificance = intChisSq;
% intStats.coeffs = intCoeffs;
stats.chiSqSignificance = chisq;
% stats.coeffs = coeffs;
end




function controlledGenoPart = controlGenoCovariates(genoInt, covariates)
controlledGenoPart = genoInt;
ppb = ParforProgressbar(length(genoInt));
parfor i=1:length(genoInt)
    controlledGenoPart{i} = single(nan * zeros(size(genoInt{i})));
    selection = ~any(genoInt{i} == 255,2);

    controlledGenoPart{i}(selection, :) = single(getResiduals(covariates(selection, :),  single(genoInt{i}(selection, :))));
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












