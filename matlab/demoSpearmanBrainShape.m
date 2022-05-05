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
outDir = ['../results/' MODALITY '/spearman_correlation/' DATASET_ID '/brain_shape/'];
    
if ~isfile([outDir 'spearman_corr.csv']) || RECOMPUTE
    %%
    load(['../results/asymmetry/meta_analysis/' DATASET_ID '/par.mat']);
    %%
    snps = readtable('../SAMPLE_DATA/w_hm3.noMHC.snplist','FileType','text');
    %%
    load('../results/correlation_sources/BRAIN_SHAPE_CHR.mat');
    %%
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
end
spearman_corr = readtable([outDir 'spearman_corr.csv']);
GC = table2array(spearman_corr);
spearman_pcorr = readtable([outDir 'spearman_pcorr.csv']);
pGC = table2array(spearman_pcorr);
spearman_secorr = readtable([outDir 'spearman_secorr.csv']);
seGC = table2array(spearman_secorr);
%%

traits = spearman_corr.Properties.VariableNames;
traitsNames = cellfun(@(x)str2double(x(4:end)), traits);
traits = cellfun(@(x)['shape' x(4:end)], traits, 'UniformOutput',false);
[~,a] = sort(traitsNames);
[fig, ax] = makeHeatmap(GC(:,a), traits(a),0, nan, 0);
pos = get(fig, 'Position');
pos(4) = pos(4)/4;
set(fig, 'Position', pos);
set(ax, 'xtick', 1:10:length(traits));
xticklabels(ax, traits(a(1:10:length(traits))));

saveas(fig, [outDir 'brainShapeHeatmap.svg'])
[fig, ax] = makeHeatmap(pGC, traits,1, get(fig, 'Position'), 0);
set(ax, 'xtick', 1:10:length(traits));
xticklabels(ax, traits(a(1:10:length(traits))));
saveas(fig, [outDir 'brainShapeHeatmapPvalues.svg'])
