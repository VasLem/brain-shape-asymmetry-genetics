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

%%
DATASET_ID= [ 'STAGE00DATA/mean_imputed/' SUBSAMPLED_ID];
load(['../results/' MODALITY '/meta_analysis/' DATASET_ID '/par.mat']);
PAR1 = PAR;
DATASET_ID=['BATCH2_2021_DATA/mean_imputed/' SUBSAMPLED_ID];
load(['../results/' MODALITY '/meta_analysis/' DATASET_ID '/par.mat']);
PAR2 = PAR;
%%
outDir = ['../results/' MODALITY '/spearman_correlation/' DATASET_ID '/discoveryVReplication/'];
if ~isfolder(outDir), mkdir(outDir); end
%%
load('../results/w_hm3.noMHC.mat')
%%
[ind12, ind21] = vlookupFast(PAR1.RS, PAR2.RS);

par1P = PAR1.P(ind21,:);
par1RS = PAR1.RS(ind21);
par2P = PAR2.P(ind12, :);
%%
[ind12,ind21] = vlookupFast(par1RS, SNP_NOMHC.RS);
par1P = par1P(ind21,:);
par2P = par2P(ind21,:);
chr = SNP_NOMHC.CHR(ind12);
pos = SNP_NOMHC.POS(ind12);
%%
test = load('../SAMPLE_DATA/ldblocks_euro.mat');
%%
test.REF.LDblocks.n = length(test.REF.LDblocks.CHR);
%%
[GC,pGC,seGC] = geneticCorrelationLDblocksWithSE(chr,pos,par1P,par2P,test.REF.LDblocks,'mean');
%%
figure;
imagesc(GC)

%%
colorbar

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
saveas(fig, [outDir 'heatmap.svg'])
