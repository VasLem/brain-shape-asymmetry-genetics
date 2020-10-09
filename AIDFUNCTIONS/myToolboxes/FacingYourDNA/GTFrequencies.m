function [fAA,faA,faa] = GTFrequencies(SNP)
         N = 2*ones(size(SNP));N(isnan(SNP)) = 0;N = nansum(N,1);
         fAA = SNP;fAA(SNP==-1) = 1;fAA(SNP==0) = 0;fAA(SNP==1)=0;fAA = nansum(fAA,1)./N;
         faA = SNP;faA(SNP==-1) = 0;faA(SNP==0) = 1;faA(SNP==1)=0;faA = nansum(faA,1)./N;
         faa = SNP;faa(SNP==-1) = 0;faa(SNP==0) = 0;faa(SNP==1)=1;faa = nansum(faa,1)./N;
end