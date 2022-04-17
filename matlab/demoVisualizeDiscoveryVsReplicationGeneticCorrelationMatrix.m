clear, close all
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
out_dir = ['../results/' MODALITY '/ldsc/disVsRep/rg/'];
inp_file = [out_dir 'discovery_vs_replication_correlation.csv'];
inp = readtable(inp_file);

if iscellstr(inp.rg)
    na_mask = cellfun(@(x)(strcmp(x,'NA')), inp.rg);
    inp.rg(na_mask) = {'nan'};
    inp.rg = cellfun(@(x)(str2double(x)), inp.rg);
end
ret = inp.rg';
clrLim = [min(min(ret)),max(max(ret))]; 
diamLim = [0.3, 1];
fig=figure();
[nr,nc] = size(ret);
ImageAlpha = ~isnan(ret);
imagesc(ret, "AlphaData",ImageAlpha)
colormap(gca,'jet');
colorbar();
caxis(clrLim);
yticks([]);
axis tight
saveas(fig, [out_dir 'dirVsRepHeatmap.svg'])
