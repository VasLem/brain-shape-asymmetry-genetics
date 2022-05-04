clear all
addpath(genpath('AIDFUNCTIONS'));
%%
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('BrainAsymmetrySignificanceAnalysis'));
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'LMU Serif');
%%
DATASET = 'joinedDatasets';
REDUCTION = 1;
MODALITY = 'asymmetry';
switch REDUCTION
    case 1
    SUBSAMPLED_ID = 'not_subsampled';
    case 10
    SUBSAMPLED_ID = 'not_subsampled';
    otherwise
        error("REDUCTION=1 or 10")
end
DATASET_ID= [DATASET '/mean_imputed/' SUBSAMPLED_ID];
%%
load(['../results/asymmetry/meta_analysis/' DATASET_ID '/par.mat']);
%%
snps = readtable('../SAMPLE_DATA/w_hm3.noMHC.snplist','FileType','text');
%%
load('../results/correlation_sources/BRAIN_SHAPE_CHR.mat');
%%
outDir = ['../results/' MODALITY '/spearman_correlation/' DATASET_ID '/brain_shape/'];
if ~isfolder(outDir), mkdir(outDir); end
%%
load('../results/w_hm3.noMHC.mat')
%%
[ind12, ind21] = vlookupFast(PAR.RS, SYN.RS);

parP = PAR.P(ind21,:);
parRS = PAR.RS(ind21);
synP = SYN.P(ind12, :);
%%
[ind12,ind21] = vlookupFast(parRS, SNP_NOMHC.RS);
parP = parP(ind21,:);
synP = synP(ind21,:);
chr = SNP_NOMHC.CHR(ind12);
pos = SNP_NOMHC.POS(ind12);
%%
test = load('../SAMPLE_DATA/ldblocks_euro.mat');
%%
test.REF.LDblocks.n = length(test.REF.LDblocks.CHR);
%%
[GC,pGC,seGC] = geneticCorrelationLDblocksWithSE(chr,pos,parP,synP,test.REF.LDblocks,'mean');
%%
spearman_corr = array2table(GC,'VariableNames',SYN.NAMES);
writetable(spearman_corr, [outDir 'spearman_corr.csv']);
spearman_pcorr = array2table(pGC,'VariableNames',SYN.NAMES);
writetable(spearman_pcorr, [outDir 'spearman_pcorr.csv']);
spearman_secorr = array2table(seGC,'VariableNames',SYN.NAMES);
writetable(spearman_secorr, [outDir 'spearman_secorr.csv']);
%% 

diamLim = [0.3, 1];
traits = SYN.NAMES;
ret = GC;
ret(pGC>0.05) = nan;
clrLim = [nanmin(ret,[],'all'),nanmax(ret,[],'all')];
fig=figure;
[nr,nc] = size(ret);
imagesc(ret, 'AlphaData',~isnan(ret))



colormap(gca,'jet');
colorbar();
caxis(clrLim);
set(gca, 'YMinorTick','on')
axis equal
axis tight
saveas(fig, [outDir 'brainShapeHeatmap.svg'])
