clear all

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
load(['../results/asymmetry/meta_analysis/' DATASET_ID '/mean_imputed/not_subsampled/par.mat']);
load('../results/w_hm3.noMHC.mat','SNP_NOMHC')
%%
load('../results/correlation_sources/OTHER_TRAITS_GWAS.mat');
%%
outDir = ['../results/' MODALITY '/spearman_correlation/' DATASET_ID '/other_traits/'];
if ~isfolder(outDir), mkdir(outDir); end
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
fig=figure(Position=[0,0,400,1000]);
[nr,nc] = size(ret);
imagesc(ret, 'AlphaData',~isnan(ret))
colormap(gca,'jet');
colorbar();
caxis(clrLim);
set(gca,'xtick',1:length(traits));
set(gca, 'YMinorTick','on')
xticklabels(traits);
axis equal
axis tight
saveas(fig, [outDir 'otherTraitsHeatmap.svg'])
%%
featMats{1} = round(100 * ret)/100;
featMatsIds{1} = 'otherTraits';
featsClassesNames = traits;
datasetName = DATASET;
drawFeaturesOnPolarPartitionsGraph(featMats, featMatsIds, featsClassesNames, MODALITY, outDir, REDUCTION, 1)

