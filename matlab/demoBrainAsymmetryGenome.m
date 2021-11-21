close all;clear;
restoredefaultpath;
addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('/opt/SNPLIB/'));
rmpath('/opt/SNPLIB/mexfiles/');% to remove the functions that are causing matlab to crash
addpath(genpath('SNPLIB-master/mexfiles/'))% where I stored the re-mexed files
%%
DATA_DIR = '../SAMPLE_DATA/';


THREADS= 12;
MAX_NUM_FEATS = 0;
DATASET_INDEX = 1;

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
disp("Loading phenotype and covariates..")
pheno = load(['../results/hierarchicalClusteringDemo/' DATASET_NAME '/asymmetry_reduction10/ccPriorSegmentation/levels4_mine/phenotype_varThres80.mat']);
covariates = load(covGenoPath).COV;
disp("Initializing genomic analysis..")
try
    parpool('local',THREADS);
catch
end
%%
ST_CHR=1;
EN_CHR=22;
%%
for CHR=ST_CHR:EN_CHR
    %%
    disp(['CHR:' , num2str(CHR)]);
    GENE_SET_METHOD = 'perSNP';
    CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
    if ~isfolder(CHR_DIR), mkdir(CHR_DIR); end
    try
        load([CHR_DIR 'plink_data.mat'])
        disp("Genomic Data loaded from disk")
    catch
        disp("Loading PLINK preprocessed data..")
        obj = SNPLIB();
        obj.nThreads = THREADS;
        [snps, samples] = obj.importPLINKDATA(['../SAMPLE_DATA/IMAGEN/BRAIN/' BIOBANK '/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_' GENO_ID '_chr' num2str(CHR)]);
        geno = obj.UnpackGeno();
        geno(geno==-1) = 255;
        geno = uint8(geno);
        disp("Computing LD scores..") 
        ld = obj.CalcLDScores(snps);
        disp("Computing Allele Frequencies..")
        af = obj.CalcAlleleFrequency();
        disp("Saving to disk..")
        clear obj;
        save([CHR_DIR 'plink_data.mat'], 'af', 'ld', 'geno', 'snps', 'samples', '-v7.3')
    end
    %%
    %% Align covariates with genotype
    disp("Aligning covariates to genotype..");
    covAssignmentMatrix = (str2double(samples.IID) == str2double(covariates.IID)');
    [covGenoIndex, covIndex] = find(covAssignmentMatrix);
    clear covAssignmentMatrix
    covData = covariates.DATA(covIndex, :);
    geno =geno(covGenoIndex, :);
    samples = samples(covGenoIndex, :);
    
    
    %%
    genoId = str2double(samples.IID);
    phenoId = str2double(pheno.preprocPhenoIID);
    
    
    %% Align phenotype with genotype
    disp("Aligning genotype to phenotype..");
    assignmentMatrix = (phenoId == genoId');
    [phenoIndex, genoIndex] = find(assignmentMatrix);
    geno = geno(genoIndex, :);
    covData = covData(genoIndex,:);
    samples = samples(genoIndex, :);
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
    %% Partition SNPs
    try
        load([CHR_DIR 'intervals.mat'], 'genoInt', 'intervals', 'gId');
        disp("Loaded intervals..");
    catch
        disp("Computing Intervals..")
        intervals = getIntervals(snpsPruned, GENE_SET_METHOD);
        disp("Splitting genome into groups, according to intervals information..")
        gId = zeros(intervals(end,2),1);
        gId(intervals(:, 1)) = 1;
        gId = cumsum(gId);
        genoInt = splitapply( @(x){x'}, genoPruned', gId );
        save([CHR_DIR 'intervals.mat'], 'genoInt', 'intervals', 'gId', '-v7.3');
    end
    clear genoPruned
    %%
    
    try
        load([CHR_DIR 'controlled_geno.mat'], 'genoControlledInt');
        disp("Loaded controlled genome for covariates based on intervals.")
    catch
         disp("Controlling genome for covariates based on intervals..")
        genoControlledInt = controlGenoCovariates(genoInt, covData);
        save([CHR_DIR 'controlled_geno.mat'], 'genoControlledInt', '-v7.3');
    end
    clear genoInt
    %%
    
    try
        load([CHR_DIR 'noPartCCA.mat'], 'noPartitionStats', 'noPartitionIntStats');
        disp("Loaded Computed CCA without phenotypic partitioning");
    catch
        disp("Computing CCA without phenotypic partitioning..")
        % Without Taking Partitioning Into Consideration
        noPartPheno = pheno.clusterPCAPhenoFeatures{1};
        if MAX_NUM_FEATS~=0
            noPartPheno = noPartPheno(:,1:MAX_NUM_FEATS);
        end
        [noPartitionStats, noPartitionIntStats] = runCCA(noPartPheno, genoControlledInt, intervals, gId);
        save([CHR_DIR 'noPartCCA.mat'], 'noPartitionStats', 'noPartitionIntStats', '-v7.3');
        clear noPartPheno;
    end
    %%
    noPartitionThres = 5*10^-8;
    plotSimpleGWAS(intervals, noPartitionIntStats, CHR, noPartitionThres,  [CHR_DIR  'noPartition_feats' num2str(MAX_NUM_FEATS)]);
    %%
    gTLPartsPThres = noPartitionThres / length(pheno.clusterPCAPhenoFeatures);
    try
        load([CHR_DIR 'withPartCCA.mat'],  'gTLPartStats', 'gTLPartIntStats');
        disp("Loaded Computed CCA with phenotypic partinioning");
    catch
        
        % With Global-To-Local Partitioning Into Consideration
        disp("Computing CCA with phenotypic partitioning..")
        [gTLPartStats, gTLPartIntStats] = runCCAOnEachPartition(pheno, genoControlledInt, intervals, gId, CHR_DIR, MAX_NUM_FEATS);
        
        plotPartitionsGWAS(intervals, gTLPartIntStats, CHR, gTLPartsPThres, [CHR_DIR 'PartitionedGTL_feats' num2str(MAX_NUM_FEATS)]);
        save([CHR_DIR, 'withPartCCA.mat'], 'gTLPartStats', 'gTLPartIntStats', '-v7.3');
    end
    
    %% Significant SNPs tables extraction
    disp("Extracting significant SNPs tables..")
    %%
    prepareSignificantTablesOnEachPartition(snpsPruned,  gTLPartIntStats , intervals, gTLPartsPThres, [CHR_DIR, 'PartitionedGTLWithBC_feats' num2str(MAX_NUM_FEATS)]);
    %%
    prepareSignificantTablesOnEachPartition(snpsPruned,  gTLPartIntStats , intervals, noPartitionThres, [CHR_DIR, 'PartitionedGTLWoutBC_feats' num2str(MAX_NUM_FEATS)]);
%     
    
    %%
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
legend;
title(['Chromosome ' num2str(chromosome) ', Partitions: ' num2str(pNum)])
ylabel('-log10p');
yline(-log10(pThres));
%%
savefig(fig, [path '_logPlot.fig']);
saveas(fig, [path '_logPlot.png']);
end


function fig = plotSimpleGWAS(intervals, intStats, chromosome, pThres, path)
%%
fig = figure;
scatter(intervals(:, 1), -log10(intStats.chiSqSignificance),'.');
title(['Chromosome ' num2str(chromosome) ', ' num2str(sum(intStats.chiSqSignificance<pThres)) ' significant intervals out of ' num2str(length(intStats.chiSqSignificance))])
ylabel('-log10p');
yline(-log10(pThres));
%%
savefig(fig, [path '_logPlot.fig']);
saveas(fig, [path '_logPlot.png']);
end

function [intSigSnps, sigSnps] = prepareSignificantTablesOnEachPartition(snps, ccaIntStats, ccaIntervals, pThres, save_path)
pNum = size(ccaIntStats.chiSqSignificance ,1);
sigSnps= [];
for i=1:pNum
    partIntStats.chiSqSignificance = ccaIntStats.chiSqSignificance(i, :);
%     partStats.coeffs = ccaStats.coeffs(i, :);
    [intRet, ret] =  prepareSignificantTables(snps, partIntStats, ccaIntervals, pThres);
    if isempty(ret), continue; end
    intRet.PARTITION = repmat(i,height(intRet),1);
    ret.PARTITION = repmat(i,height(ret),1);
    if isempty(sigSnps)
        intSigSnps = intRet;
        sigSnps = ret;
    else
        if size(ret,1) ~= 0
            intSigSnps = [intSigSnps;intRet];
            try
                sigSnps = [sigSnps;ret];
            catch
                disp('h');
            end
        end
    end
end
if ~isempty(sigSnps)
    writetable(sigSnps, [save_path 'significant_snps.csv']);
    writetable(intSigSnps, [save_path 'significant_intervals.csv']);
end
end


function [intSigSnps, sigSnps] = prepareSignificantTables(snps, ccaIntStats, ccaIntervals, pThres, savePath)
arguments
    snps
    ccaIntStats
    ccaIntervals
    pThres
    savePath string = ""
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

function [stats, intStats, intervals] = runCCAOnEachPartition(pheno, geno, intervals, intIdVec, chr_dir, maxNumPhenoFeats)
arguments
    pheno
    geno
    intervals
    intIdVec
    chr_dir
    maxNumPhenoFeats =0;
    
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
    PART_DIR = [chr_dir 'par' num2str(k) '/'];
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
arguments
    genoInt cell
    covariates double
end
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
arguments
    snps table
    geneSetMethod string='perSNP'
    namedArgs.windowSize double=20000
end
switch geneSetMethod
    case 'slidingWindow'
        snpCnt = 1;
        cnt = 1;
        intSlidWindow = zeros(height(snps),2);
        while snpCnt<height(snps)
            startInd = snpCnt;
            startPos = snps.POS(snpCnt);
            offset = find(snps.POS(snpCnt:end) > startPos + namedArgs.windowSize, 1,  'first');
            if size(offset)==0, offset = height(snps) + 1; end
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
        intervals(:, 2) = [ia(2:end)-1;height(snps.POS)];
end
end












