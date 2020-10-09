function OUT = myReduceGWAS(GWAS,index)
         if isfield(GWAS,'CHR'), OUT.CHR = GWAS.CHR(index); end
         if isfield(GWAS,'POS'), OUT.POS = GWAS.POS(index); end
         if isfield(GWAS,'RS'), OUT.RS = GWAS.RS(index); end
         if isfield(GWAS,'JOBID'), OUT.JOBID = GWAS.JOBID(index); end
         if isfield(GWAS,'A1'), OUT.A1 = GWAS.A1(index); end
         if isfield(GWAS,'A2'), OUT.A2 = GWAS.A2(index); end
         if isfield(GWAS,'MAF'), OUT.MAF = GWAS.MAF(index); end
         if isfield(GWAS,'N'), OUT.N = GWAS.N(index,:); end
         if isfield(GWAS,'PLH'), OUT.PLH = GWAS.PLH(index,:); end
         if isfield(GWAS,'PRH'), OUT.PRH = GWAS.PRH(index,:); end
         if isfield(GWAS,'P')
             switch ndims(GWAS.P)
                 case 2
                     OUT.P = GWAS.P(index,:); 
                 case 3
                     OUT.P = GWAS.P(index,:,:); 
             end
         end
         if isfield(GWAS,'SNP'), OUT.SNP = GWAS.SNP(index,:); end
         if isfield(GWAS,'bestP'), OUT.bestP = GWAS.bestP(index); end
         if isfield(GWAS,'Study'), OUT.Study = GWAS.Study(index); end
         if isfield(GWAS,'PSH'), OUT.PSH = GWAS.PSH(index,:); end
         if isfield(GWAS,'PSWH'), OUT.PSWH = GWAS.PSWH(index,:); end
         if isfield(GWAS,'CHISH'), OUT.CHISH = GWAS.CHISH(index,:); end
         if isfield(GWAS,'CHISWH'), OUT.CHISWH = GWAS.CHISWH(index,:); end
         if isfield(GWAS,'CHI'), OUT.CHI = GWAS.CHI(index,:); end
end