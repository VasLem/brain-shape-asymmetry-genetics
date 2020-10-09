function out = getMAF(geno)
         geno = codetoM(geno);
         tmp = geno(geno>=0);
         N = length(tmp)*2;
         tmp = 2-tmp;
         Na = sum(tmp);
         out = Na./N;
end