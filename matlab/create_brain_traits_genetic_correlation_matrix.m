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
labels = {'auto', 'face', 'brain_shape'};
for trait_ind=1:3
    switch trait_ind
        case 1
            trait='';
        case 2
            trait='face';
        case 3
            trait='brain_shape';
    end
    out_dir = ['../results/ldsc/' DATASET_ID '/rg/'];
    if trait_ind==1
        inp_file = [out_dir '/correlation.csv'];
    else
        inp_file = [out_dir '/' trait '_correlation.csv'];
    end
    inp = readtable(inp_file);

    if trait_ind==1
        ret = eye(31);
        ret(ret==0) = nan;
    else
        ret = nan(31,round(height(inp)/31));
    end
    if iscellstr(inp.rg)
        na_mask = cellfun(@(x)(strcmp(x,'NA')), inp.rg);
        inp.rg(na_mask) = {'nan'};
        inp.rg = cellfun(@(x)(str2double(x)), inp.rg);
        na_mask = cellfun(@(x)(strcmp(x,'NA')), inp.p);
        inp.p(na_mask) = {'nan'};
        inp.p = cellfun(@(x)(str2double(x)), inp.p);
    end
    for row=1:height(inp)
        if isnan(inp.p) | (inp.p(row)>0.05)
            continue
        end
        [~,fname1] = fileparts(inp.p1(row));
        [~,fname2] = fileparts(inp.p2(row));

        ind1 = str2double(regexprep(fname1, "par(\d+).sumstats", '$1'));
        ind2 = str2double(regexprep(fname2, "par(\d+).sumstats", '$1'));
        ret(ind1, ind2) = inp.rg(row);
        if trait_ind==1
            ret(ind2, ind1) = inp.rg(row);
        end
    end
    clrLim = [min(min(ret)),max(max(ret))]; 
    diamLim = [0.3, 1];
    fig=figure();
    imagesc(ret)
    colormap(gca,'jet');
    colorbar();
    caxis(clrLim);
    axis equal
    axis tight
    saveas(fig, [out_dir labels{trait_ind} '_heatmap.svg'])
end