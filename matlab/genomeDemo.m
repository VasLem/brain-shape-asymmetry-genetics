addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('/opt/SNPLIB/'));
rmpath('/opt/SNPLIB/mexfiles/');% to remove the functions that are causing matlab to crash
addpath(genpath('SNPLIB-master/mexfiles/'))% where I stored the re-mexed files
chr = 17;
try
    load(['../results/genomics/data_chr' num2str(chr) '.mat'])
catch
    %%
    obj = SNPLIB();
    obj.nThreads = 16;
    [snps, samples] = obj.importPLINKDATA(['../SAMPLE_DATA/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_sel19908_chr' num2str(chr)]);
    geno = uint8(obj.UnpackGeno());
    af = obj.CalcAlleleFrequency();
    save(['../results/genomics/data_chr' num2str(chr) '.mat'], 'af', 'geno', 'snps', 'samples', '-v7.3')
end
%%
RESULTS_DIR = '../results/genomeDemo/';
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
%%
genoId = str2double(samples.('IID'));

%%
pheno = load('../results/hierarchicalClusteringDemo/AsymmetryPhenotype.mat');
phenoId = str2double(pheno.phenoIID);


%% Allign phenotype with genotype
assignmentMatrix = (phenoId == genoId');
[phenoIndex, genoIndex] = find(assignmentMatrix);
geno = geno(genoIndex, :);


%% Represent SNPs that are not really SNPs (one replaced by many, or many replaced by one)
[sortedAf, sortedAfArg] = sort(af);
sortedSnps = snps(sortedAfArg, :);
%% Decide whether I should remove multi-character allele correspondences (Incorporating Indels)
close all
fig = figure("Name", "Occurences of Multi-Character Alleles");
hold on
yyaxis left
lHandle = plot(sortedAf);
set(lHandle, {'DisplayName'}, {'Ordered Frequency'});

oneToMany = strlength(sortedSnps.(4))>1 & strlength(sortedSnps.(5)) == 1;
manyToMany = strlength(sortedSnps.(4))>1 & strlength(sortedSnps.(5)) > 1;
manyToOne = strlength(sortedSnps.(4))==1 & strlength(sortedSnps.(5)) > 1;
edges = linspace(0, size(manyToOne ,1), 50);
oTMH = histcounts(find(oneToMany),edges);
mTOH = histcounts(find(manyToOne),edges);


sHandle = scatter(find(manyToMany), 0.95 * ones(1, sum(manyToMany))', 'x');
sHandle.DisplayName = ['N->N: ', num2str(sum(manyToMany))];

yyaxis right
bHandle = bar( 0.5*(edges(2:end) + edges(1:end-1)), [oTMH', mTOH'], 1);
bHandle(1).FaceColor = 'blue';
bHandle(1).DisplayName = ['1->N:  '  num2str(sum(oneToMany))];
bHandle(1).EdgeColor = 'k';
bHandle(2).FaceColor = 'red';ax = gca;
bHandle(2).DisplayName = ['N->1 :', num2str(sum(manyToOne))];
bHandle(2).EdgeColor = 'k';
xlabel("Sorted Snps")
ax = gca;
ax.YAxis(2).Color = 'k';

legend;
hold off
savefig(fig, [RESULTS_DIR 'multi_characters_alleles.fig']);
saveas(fig, [RESULTS_DIR 'multi_characters_alleles.png']);


%% Remove multiallelic for now.
indelFilter = strlength(snps.(4))==1 & strlength(snps.(5))==1;
genoPruned = geno(:, indelFilter);
snpsPruned = snps(indelFilter, :);


%% Identify the genetic indices for genes with more than one allele
[rsids , ~, ic] = unique(snpsPruned.('RSID'));
alleleCounts = accumarray(ic,1);
%% Some checks
%regarding the fact that positions need to be sorted
assert(all(sort(snpsPruned.POS) == snpsPruned.POS));
% and that all phenotype ids  correspond to genotype ones
assert(all((1:length(phenoId)) == phenoIndex'));
%% Partition SNPs
% Filter candidate SNPs first by testing them against the whole
% phenotype, using means
%
% meanFeats = pheno.means;
% canoncorr(meanFeats);

% Across all brain segments
% phenoPart = pheno.clusterPCAPhenoFeatures{clusterInd};
phenoWhole = horzcat(pheno.clusterPCAPhenoFeatures{:});

[stats, intStats, intervals] = runCCA(phenoWhole, genoPruned, snpsPruned);


%%
pThres = 0.05/ size(intStats.chiSqSignificance, 1);
significantInts = intervals(intStats.chiSqSignificance<pThres, :);
starts = significantInts(:,1);
ends = significantInts(:,2);
stSnps = snpsPruned(starts,:);
stSnps.Properties.VariableNames = strcat(snpsPruned.Properties.VariableNames', '_ST')';
enSnps = snpsPruned(ends,:);
enSnps.Properties.VariableNames = strcat(snpsPruned.Properties.VariableNames', '_EN')';

sigSnpsCoef = zeros(size(significantInts, 1),1);
sigSnps = table();
for c=1:size(significantInts,1)
    intSnps = snpsPruned(starts(c):ends(c),:);
    intCoeffs = stats.coeffs(starts(c):ends(c));
    [~, maxCoefArg] = max(abs(intCoeffs));
    sigSnpsCoef(c) = intCoeffs(maxCoefArg);
    sigSnps(c, :) = intSnps(maxCoefArg, :);
end
sigSnps.COEF = sigSnpsCoef;

writetable(sigSnps, [RESULTS_DIR, 'significant_snps_chr' num2str(chr) '_slidwin_20k.csv']);
intTable = [stSnps, enSnps];
writetable(intTable, [RESULTS_DIR, 'significant_intervals_chr' num2str(chr) '_slidwin_20k.csv']);


%%
fig = figure;
scatter(intervals(:, 1), -log10(intStats.chiSqSignificance),'.');
title(['Chromosome ' num2str(chr) ', ' num2str(sum(intStats.chiSqSignificance<pThres)) ' significant intervals out of ' num2str(length(intStats.chiSqSignificance))])
ylabel('-log10p');
yline(-log10(pThres));
%%
savefig(fig, [RESULTS_DIR 'chr' num2str(chr) '_all_segments.fig']);
saveas(fig, [RESULTS_DIR 'chr' num2str(chr) '_all_segments.png']);


%%

function [stats, intStats, intervals] = runCCA(phenoPart, geno, snps, method, namedArgs)
arguments
    phenoPart double
    geno uint8
    snps table
    method string='slidingWindow'
    namedArgs.windowSize double=20000
end

snpCnt = 1;
cnt = 1;
if strcmp(method, 'slidingWindow')
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
    intervals = intSlidWindow;
end
% haplotype based

%%

intChiSqSignificance = zeros(size(intervals, 1),1);
[path,ID] = setupParForProgress(length(intervals));
parfor i=1:length(intervals)
    startInd = intervals(i, 1);
    endInd = intervals(i, 2);
    
    genoPart = single(geno(:, startInd:endInd));
    [A,~,~,~,~,st] = canoncorr(genoPart, phenoPart);
    intChiSqSignificance(i) = st.pChisq(1);
    intCoeffs{i} = A(:,1);
    parfor_progress;
end
closeParForProgress(path,ID);
intStats.chiSqSignificance = intChiSqSignificance;
intStats.coeffs = intCoeffs;

%%
chiSqSignificance = zeros(size(geno,2),1);
coeffs = zeros(size(geno,2), 1);
for i=1:length(intervals)
    startInd = intervals(i, 1);
    endInd = intervals(i, 2);
    chiSqSignificance(startInd:endInd) = intChiSqSignificance(i);
    coeffs(startInd:endInd) = intCoeffs{i};
end
stats.chiSqSignificance = chiSqSignificance;
stats.coeffs = coeffs;
end














