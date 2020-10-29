function [fA,fa] = getAlleleFrequencies(SNP)
         [s,~] = size(SNP);
         N = 2*ones(size(SNP));N(isnan(SNP)) = 0;N = nansum(N,1);
         fA = SNP;fA(SNP==-1) = 2;fA(SNP==0) = 1;fA(SNP==1) = 0;fA = nansum(fA,1)./N;
         fa = SNP;fa(SNP==-1) = 0;fa(SNP==0) = 1;fa(SNP==1) = 2;fa = nansum(fa,1)./N;
end