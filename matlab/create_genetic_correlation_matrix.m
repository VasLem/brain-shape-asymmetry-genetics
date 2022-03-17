DATASET='STAGE00DATA';
for trait_ind=1:3
    switch train_ind
        case 1
            trait='';
        case 2
            trait='face';
        case 3
            trait='brain_shape';
    end
    out_dir = ['../results/ldsc/' DATASET ];
    inp_file = [out_dir '/' trait 'correlation.csv'];
    inp = readtable(inp_file);
    %%
    if trait_ind==1
        ret = eye(31);
    else
        ret = zeros(31,height(inp)/31);
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
    clrLim = [0,1]; 
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