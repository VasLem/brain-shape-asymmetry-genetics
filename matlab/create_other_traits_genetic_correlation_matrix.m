clear, close all
DATASET = 'STAGE00DATA';
REDUCTION = 1;
switch REDUCTION
    case 1
    SUBSAMPLED_ID = 'not_subsampled';
    case 10
    SUBSAMPLED_ID = 'not_subsampled';
    otherwise
        error("REDUCTION=1 or 10")
end
DATASET_ID= [DATASET '/mean_imputed/' SUBSAMPLED_ID];
out_dir = ['../results/ldsc/' DATASET_ID '/rg/'];
inp_file = [out_dir '/other_traits_correlation.csv'];
inp = readtable(inp_file);

ret = nan(31,round(height(inp)/31));
if iscellstr(inp.rg)
    na_mask = cellfun(@(x)(strcmp(x,'NA')), inp.rg);
    inp.rg(na_mask) = {'nan'};
    inp.rg = cellfun(@(x)(str2double(x)), inp.rg);
end
traits = ["ADHD" ,"ASD", "BipolarDisorderAll", "Handedness", "Intelligence", "LanguageFunctionalConnectivity",  "OCD",  "Schizophrenia"];
for row=1:height(inp)
    if inp.p(row) > 0.05
        continue
    end
    [~,fname1] = fileparts(inp.p1(row));
    fname2 = inp.p2(row);

    ind1 = str2double(regexprep(fname1, "par(\d+).sumstats", '$1'));
    trait =  regexprep(fname2, "../SAMPLE_DATA/OTHER_TRAITS_GWAS/(\w+)/munged.sumstats.gz", '$1');
    ind2 =  find(contains(traits,trait));
    ret(ind1, ind2) = inp.rg(row);
end
clrLim = [min(min(ret)),max(max(ret))]; 
diamLim = [0.3, 1];
fig=figure();
imagesc(ret)
colormap(gca,'jet');
colorbar();
caxis(clrLim);
xticklabels(traits);
axis equal
axis tight
saveas(fig, [out_dir 'otherTraitsHeatmap.svg'])
featMats{1} = round(100 * ret)/100;
featMatsIds{1} = 'otherTraits';
featsClassesNames = traits;
datasetName = DATASET;
drawFeaturesOnPolarPartitionsGraph(featMats, featMatsIds, featsClassesNames, datasetName, out_dir, REDUCTION, 0)