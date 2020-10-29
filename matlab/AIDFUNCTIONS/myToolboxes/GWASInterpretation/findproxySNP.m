function out = findproxySNP(chr,pos,GWAS,distnet,mindist,ldT)
         n = length(PEAKS.POS);
         exactOverlap = zeros(n,1);
         proxOverlap = zeros(n,1);
         for i=1:1:n
            %i=1;
            disp(num2str(i));
            index = find(GWAS.CHR==PEAKS.CHR(i));
            OUT = myReduceGWAS(GWAS,index); 
            index = find(uint32(OUT.POS)==uint32(PEAKS.POS(i)));
            if ~isempty(index)
               exactOverlap(i) = 1;
            end
            index = find(OUT.POS<=(PEAKS.POS(i)+distnet)&OUT.POS>=(PEAKS.POS(i)-distnet));
            if ~isempty(index)
                OUT = myReduceGWAS(OUT,index);
                GENO2 = PEAKS.SNP(i,:);
                GENO1 = OUT.SNP;
                R2 = getLDStructure(GENO1,GENO2);
                R2(isnan(R2)) = 0;
                index1 = find(R2>=ldT);
                posdist = pdist2(double(PEAKS.POS(i)),double(OUT.POS));
                index2 = find(posdist<=mindist);
                index = union(index1,index2);
                if ~isempty(index)
                    proxOverlap(i) = 1;
                end
            end
         end
         out.exactOverlap = exactOverlap;
         out.proxOverlap =proxOverlap;
end

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
         if isfield(GWAS,'P'), OUT.P = GWAS.P(index,:); end
         if isfield(GWAS,'SNP'), OUT.SNP = GWAS.SNP(index,:); end
         if isfield(GWAS,'bestP'), OUT.bestP = GWAS.bestP(index); end
         if isfield(GWAS,'Study'), OUT.Study = GWAS.Study(index); end
end






