%%read brain shape GWAS
load("../SAMPLE_DATA/BRAIN_SHAPE_GWAS/BRAINGWAS.mat");
%%
OUTPUT_DIR = '../SAMPLE_DATA/BRAIN_SHAPE_PARTITIONS/';
if ~isfolder(OUTPUT_DIR), mkdir(OUTPUT_DIR); end
for par=1:285
    fil = sprintf("%spar%02d.csv",OUTPUT_DIR, par);
    tab = table;
    tab.rsID = GWAS.RS;
    tab.N = GWAS.N;
    tab.A2 = GWAS.A2;
    tab.A1 = GWAS.A1;
    tab.P_value = 10.^-(single(GWAS.P(:,par))/10000);
    tab.P_value(tab.P_value==0) = min(tab.P_value(tab.P_value ~= 0));
    tab.ChiScore = GWAS.CHI(:,par);
    tab.chromosome = GWAS.CHR;
    writetable(tab, fil);
    gzip(fil);
    delete(fil);
end