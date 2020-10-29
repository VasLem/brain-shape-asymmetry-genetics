function [G,NSNPS] = genetic_distance(SNP)
N = size(SNP,1); % Number of individuals
nrSNPs = size(SNP,2); % Number of SNPs
SNP_Z = zeros(N,nrSNPs);
for i=1:nrSNPs
    GT = SNP(:,i);
    ind = find(~isnan(GT));
    GT_Z = zscore(GT(ind));
    SNP_Z(ind,i) = GT_Z;    
end
SNP_T = 1*(SNP_Z~=0);
NSNPS = SNP_T*SNP_T'; % Numbers of Individuals used for estimation
G = SNP_Z*SNP_Z'./NSNPS; % Genetic Relationship Matrix

end