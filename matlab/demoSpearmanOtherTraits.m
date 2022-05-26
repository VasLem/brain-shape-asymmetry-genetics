clear all
close all
%%
addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('BrainAsymmetrySignificanceAnalysis'));
setuplatex;
%%
DATASET = 'joinedDatasets';
REDUCTION = 1;
MODALITY = 'asymmetry';
RECOMPUTE = 0;
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

%%
outDir = ['../results/' MODALITY '/spearman_correlation/' DATASET_ID '/other_traits/'];
if ~isfile([outDir, 'spearman_corr.csv']) || RECOMPUTE
load(['../results/asymmetry/meta_analysis/' DATASET_ID '/par.mat']);
load('../results/w_hm3.noMHC.mat','SNP_NOMHC')
load('../results/correlation_sources/OTHER_TRAITS_GWAS.mat');
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
end
%%

spearman_corr = readtable([outDir 'spearman_corr.csv']);
GC = table2array(spearman_corr);
spearman_pcorr = readtable([outDir 'spearman_pcorr.csv']);
pGC = table2array(spearman_pcorr);
spearman_secorr = readtable([outDir 'spearman_secorr.csv']);
seGC = table2array(spearman_secorr);

traits = spearman_corr.Properties.VariableNames;

writetable(spearman_corr, [outDir 'spearman_corr.csv']);
diamLim = [0.3, 1];
ret = GC;
mask  =~all(isnan(ret),1);
if any(~mask)
    writelines(traits(~mask), [outDir 'spearman_failed']);
end
ret = ret(:,mask);
pRet = pGC;
pRet = pGC(:, mask);
traits = traits(mask);
%%
close all
fig = makeHeatmap(ret, traits, 0);
colormap('parula')
saveas(fig, [outDir 'otherTraitsHeatmap.svg'])

fig = makeHeatmap(pRet, traits,1, get(fig,'Position'));
colormap(flipud(parula))
saveas(fig, [outDir 'otherTraitsHeatmapPvalues.svg'])

%%
ret = GC;
ret(pGC>0.05) = nan;
ret = round(100 * ret)/100;
ret(ret==0) = nan;
featMats{1} = ret;
featMatsIds{1} = 'otherTraits';
featsClassesNames =spearman_corr.Properties.VariableNames;
datasetName = DATASET;
drawFeaturesOnPolarPartitionsGraph(featMats, featMatsIds, featsClassesNames, MODALITY, outDir, REDUCTION, 1)

