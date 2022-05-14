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
outDir = ['../results/' MODALITY '/spearman_correlation/' DATASET_ID '/filteredSnps/'];
%%
if ~isfolder(outDir), mkdir(outDir); end
%%
load(['../results/asymmetry/meta_analysis/' DATASET_ID '/par.mat']);
load('../results/w_hm3.noMHC.mat','SNP_NOMHC')
%%
[ind, ~] = vlookupFast(SNP_NOMHC.RS, PAR.RS);
filteredP = PAR.P(ind,:);
filteredCHR = PAR.CHR(ind);
filteredPOS = PAR.POS(ind);
filteredRS = PAR.RS(ind);
%%
removedInd = setdiff(1:size(PAR.P,1), ind);
%%
removedP = PAR.P(removedInd,:);
removedCHR = PAR.CHR(removedInd);
removedRS = PAR.RS(removedInd);
removedPOS = PAR.POS(removedInd);
%%
removedSignificantSnps = cell(1,31);
filteredSignificantSnps = cell(1,31);
for partition=1:31
    msk = removedP(:,partition) < 5e-8;
    removedSignificantSnps{partition} = table(removedP(msk,partition), removedRS(msk), removedPOS(msk), removedCHR(msk), 'VariableNames', {'P','RS','POS','CHR'});
    msk = filteredP(:,partition) < 5e-8;
    filteredSignificantSnps{partition} = table(filteredP(msk,partition), filteredRS(msk), filteredPOS(msk), filteredCHR(msk), 'VariableNames', {'P','RS','POS','CHR'});
end
%%
stats = cell(31,1);
for partition=1:31
    r = removedSignificantSnps{partition};
    f = filteredSignificantSnps{partition};
    stats{partition} = struct;
    for CHR=1:22
        stats{partition}.CHR(CHR) = CHR;
        stats{partition}.PART(CHR) = partition;
        cr = r(r.CHR == CHR,:);
        cf = f(f.CHR == CHR,:);
        stats{partition}.removedNum(CHR) = height(cr);
        stats{partition}.filteredNum(CHR) = height(cf);
        stats{partition}.removedBonfNum(CHR) = height(cr(cr.P<5e-8/31,:));
        stats{partition}.filteredBonfNum(CHR) = height(cf(cf.P<5e-8/31,:));
        
        
        if isempty(cr)
            stats{partition}.removedPMin(CHR)= nan;
            stats{partition}.removedPMean(CHR) = nan;
        else
            stats{partition}.removedPMean(CHR) = mean(cr.P);
            stats{partition}.removedPMin(CHR)= min(cr.P);
        end
        if isempty(cf)
            stats{partition}.filteredPMean(CHR) = nan;
            stats{partition}.filteredPMin(CHR)= nan;
        else
            stats{partition}.filteredPMean(CHR) = mean(cf.P);
            stats{partition}.filteredPMin(CHR)= min(cf.P);
        end
        
        
    end
end
%%
ret = cellfun(@(y)struct2table(structfun(@(x)x', y, 'UniformOutput',false)),stats, 'UniformOutput', false);
ret = cellfun(@(x)x(~isnan(x.filteredPMean) & ~isnan(x.removedPMean),:), ret, 'UniformOutput',false);
%%
removedAllSig = cellfun(@(x)x((x.removedNum ~= 0) & (x.filteredNum ==0), :), ret, 'UniformOutput',false);
removedAllSigPartitions = find(cellfun(@(x)~isempty(x), removedAllSig));
removedAllSig = removedAllSig(removedAllSigPartitions);
% removedAllSig = cat(removedAllSig{:});
removedAllBonfSig = cellfun(@(x)x((x.removedBonfNum ~= 0) & (x.filteredBonfNum ==0), :), ret, 'UniformOutput',false);
removedAllBonfSigPartitions = find(cellfun(@(x)~isempty(x), removedAllBonfSig));
removedAllBonfSig = removedAllBonfSig(removedAllBonfSigPartitions);
removedAllBonfSig = cat(1, removedAllBonfSig{:});

%%
f=figure('Position',[1000 222 1490 1116]);
hold on;
t = 0;
xticks_pos = [];
AVAILABLE_CHRS = [];
for CHR=1:22
    msk = (removedCHR == CHR);
    rPos= t + removedPOS(msk);
    scatter(rPos, -log10(removedP(msk,1)),15,[0.7,0.7,0.7],'.');
    msk = (filteredCHR == CHR);
    fPos = t + filteredPOS(msk);
    scatter(fPos, -log10(filteredP(msk,1)),15,'.');
    t = max(max(fPos), max(rPos));
    xticks_pos(end + 1) = median([fPos(:); rPos(:)]);
    AVAILABLE_CHRS(end+1) = CHR;
end
xticks(xticks_pos)
xticklabels(arrayfun(@num2str, AVAILABLE_CHRS , 'UniformOutput', 0))
yline(-log10(5e-8));
yline(-log10(5e-8/31),'--');
ylabel('-log10p');
set(gca,'TickDir','out');
saveas(f, [outDir,'filteredGwas.svg'])
%%
writetable(removedAllBonfSig, [outDir, 'removedAllBonf.csv'])
writetable(ret{1}, [outDir, 'entireStats.csv'])

