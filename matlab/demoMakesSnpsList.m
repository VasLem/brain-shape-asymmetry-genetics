addpath(genpath('AIDFUNCTIONS'));
load(['../results/asymmetry/meta_analysis/joinedDatasets/mean_imputed/not_subsampled/par.mat']);
snps = readtable('../SAMPLE_DATA/w_hm3.noMHC.snplist','FileType','text');
[ind12,ind21] = vlookupFast(PAR.RS, snps.SNP); 
SNP_NOMHC.POS = PAR.POS(ind21);
SNP_NOMHC.CHR = PAR.CHR(ind21);
SNP_NOMHC.RS = PAR.RS(ind21);
SNP_NOMHC.A1 = snps.A1(ind12);
SNP_NOMHC.A2 = snps.A2(ind12);
save('../results/w_hm3.noMHC.mat','SNP_NOMHC')
