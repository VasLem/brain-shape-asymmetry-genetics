close all
DATASET_ID='STAGE00DATA/mean_imputed/not_subsampled/rg';
for trait_ind=1:3
    switch trait_ind
        case 1
            trait='';
        case 2
            trait='face';
        case 3
            trait='brain_shape';
    end
    out_dir = ['../results/ldsc/' DATASET_ID ];
    if trait_ind==1
        inp_file = [out_dir '/correlation.csv'];
    else
        inp_file = [out_dir '/' trait '_correlation.csv'];
    end
    inp = readtable(inp_file);

    if trait_ind==1
        ret = eye(31);
    else
        ret = zeros(31,round(height(inp)/31));
    end
    if iscellstr(inp.rg)
        na_mask = cellfun(@(x)(strcmp(x,'NA')), inp.rg);
        inp.rg(na_mask) = {'nan'};
        inp.rg = cellfun(@(x)(str2double(x)), inp.rg);
    end
    for row=1:height(inp)
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
    saveas(fig, [out_dir 'trait.svg'])
end