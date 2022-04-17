% add to SHA et all GWAS a SIGN column filled with 1

tab = readtable('../SAMPLE_DATA/Asymmetry_SHA_GWAS/input_raw');
tab.SIGN = ones(height(tab),1);
tab.N = 32256 + zeros(height(tab),1);
writetable(tab, '../SAMPLE_DATA/Asymmetry_SHA_GWAS/input');