clear all
addpath(genpath('AIDFUNCTIONS'));
datasets = {'STAGE00DATA', 'BATCH2_2021_DATA','joinedDatasets'};
for datasetInd=1:3
    dataset = datasets{datasetInd};
    INPUT_DIR = ['../results/asymmetry/meta_analysis/' dataset '/mean_imputed/not_subsampled/'];
    PAR = struct;
    for part=1:31
        disp(['Partition:', num2str(part)]);
        inpStatsFile = [INPUT_DIR 'CCAPart' sprintf('%02d', part), '.csv.gz'];
        unzippedFile = gunzip(inpStatsFile);
        x = readtable(unzippedFile{1},FileType="text",Delimiter=',');
        delete(unzippedFile{1});
        if part == 1
            PAR.RS = x.rsID;
            try
                PAR.POS = x.position;
            catch
            end
            PAR.CHR = x.chromosome;
            PAR.P = x.P_value;
        else
            [~,ind21] = vlookupFast(x.rsID, PAR.RS);
            PAR.P = cat( 2, PAR.P, x.P_value(ind21));
        end
    end
    save([INPUT_DIR, 'par.mat'], 'PAR', '-v7.3');
end
%%
load('../results/asymmetry/meta_analysis/joinedDatasets/mean_imputed/not_subsampled/par.mat');
snps = readtable('../SAMPLE_DATA/w_hm3.noMHC.snplist','FileType','text');
%%
[ind12,ind21] = vlookupFast(PAR.RS, snps.SNP); 
SNP_NOMHC.POS = PAR.POS(ind21);
SNP_NOMHC.CHR = PAR.CHR(ind21);
SNP_NOMHC.RS = PAR.RS(ind21);
SNP_NOMHC.A1 = snps.A1(ind12);
SNP_NOMHC.A2 = snps.A2(ind12);
%%
save('../results/w_hm3.noMHC.mat','SNP_NOMHC')
%%
load('../results/correlation_sources/OTHER_TRAITS_GWAS.mat');
%%
[ind12, ind21] = vlookupFast(PAR.RS, SYN.RS);

parP = PAR.P(ind21,:);
parRS = PAR.RS(ind21);
synP = SYN.P(ind12, :);
%%
[ind12,ind21] = vlookupFast(parRS, SNP_NOMHC.RS);
parP = parP(ind21,:);
synP = synP(ind21,:);
chr = SNP_NOMHC.CHR(ind12);
pos = SNP_NOMHC.POS(ind12);
%%
test = load('../SAMPLE_DATA/ldblocks_euro.mat');
%%
test.REF.LDblocks.n = length(test.REF.LDblocks.CHR);
%%
[GC,pGC,seGC] = geneticCorrelationLDblocksWithSE(chr,pos,parP,synP,test.REF.LDblocks,'mean');
%%
spearman_corr = array2table(GC,'VariableNames',SYN.NAMES);
writetable(spearman_corr, '../results/asymmetry/spearman_corr.csv')
spearman_pcorr = array2table(pGC,'VariableNames',SYN.NAMES);
writetable(spearman_pcorr, '../results/asymmetry/spearman_pcorr.csv')
spearman_secorr = array2table(seGC,'VariableNames',SYN.NAMES);
writetable(spearman_secorr, '../results/asymmetry/spearman_secorr.csv')

