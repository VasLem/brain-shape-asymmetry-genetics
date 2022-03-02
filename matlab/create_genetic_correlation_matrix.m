DATASET='STAGE00DATA';
out_dir = ['../results/ldsc/' DATASET ];
inp_file = [out_dir '/correlation.csv'];
inp = readtable(inp_file);
%%
ret = eye(31);
for row=1:height(inp)
[~,fname1] = fileparts(inp.p1(row));
[~,fname2] = fileparts(inp.p2(row));

ind1 = str2double(regexprep(fname1, "par(\d+).sumstats", '$1'));
ind2 = str2double(regexprep(fname2, "par(\d+).sumstats", '$1'));
ret(ind1, ind2) = inp.rg(row);
ret(ind2, ind1) = inp.rg(row);
end
%%

% Set [min,max] value of C to scale colors
% This must span the range of your data!
clrLim = [0,1];  % or [-1,1] as in your question
% Set the  [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.3, 1];
% Plot the data using imagesc() for later comparison
fig=figure();
imagesc(ret)
colormap(gca,'jet');
colorbar();
caxis(clrLim);
axis equal
axis tight
saveas(fig, out_dir)