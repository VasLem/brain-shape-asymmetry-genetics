function out = ExtractPeakByLocationAndLDAndPCorr(SEL,pCrit,distnet,distmerge,ldT,cT)
%          if nargin<4, ldT=0; end% no LD
%          if nargin<3, distnet = 500e3;end
%          if nargin<2, pCrit = 5e-8;end
         if pCrit<1, pCrit = -log10(pCrit);end
         rs = SEL.RS;
         chr = SEL.CHR;
         pos = SEL.POS;
         p = SEL.P;
         bestp = SEL.bestP;
         SNP = SEL.SNP;
         % initiate
         nHits = length(SEL.POS);
         list = 1:nHits;counter = 0;lab = 0;
         labels = zeros(1,nHits);
         top = zeros(1,nHits);
         while ~isempty(list)&&lab<=nHits
             counter = counter +1;
             disp(num2str(counter));
             lab = lab+1;
             nlist = length(list);
             [m,mind] = max(bestp(list));
             bestind = list(mind);
             top(bestind) = 1;
             labels(bestind) = lab;
             list = setdiff(list,bestind);% take the best out;
             bestchr = chr(bestind);
             subindex = find(chr(list)==bestchr);% select snps on the same chromosome
             if isempty(subindex), continue; end % there are no more snps on the same chromosome
             subind = list(subindex);
             bestpos = pos(bestind);
             subpos = pos(subind);
             take = subpos<=(bestpos+distnet)&subpos>=(bestpos-distnet);% select snps in the same proximity
             if sum(take)==0, continue; end % there are no snps close by
             subind = subind(find(take));
             % CRITERIA 1: COMPUTING THE SPEARMAN CORRELATIONS BETWEEN SEGMENT P-VALUES
%              p1 = p(bestind,:);p2 = p(subind,:);% select snps correlated phenotypic effect
%              pcorr = zeros(1,length(subind));
%              for k=1:length(subind)
%                 pcorr(k) = corr(p1',p2(k,:)','type','spearman'); 
%              end
%              subindSC = subind(find(pcorr>=cT));
             % CRITERIA 2: LD
             GENO1 = SNP(bestind,:);
             GENO2 = SNP(subind,:);
             R2 = getLDStructure(GENO2,GENO1);
             R2(isnan(R2)) = 1;
             subindLD = subind(find(R2>=ldT));
             % CRITERIA 3: DIST
             posdist = pdist2(double(pos(bestind)),double(pos(subind)));
             subindDIST = subind(find(posdist<=distmerge));
             % MERGING CRITERIA
             subind =union(subindDIST,subindLD);
%              subind =union(subind,subindSC);
             labels(subind) = lab;
             list = setdiff(list,subind);
         end
         peakindex = find(top);
         subindex = find(bestp(peakindex,:)>=pCrit);
         peakindex = peakindex(subindex); % only retain peaks with pCRIT top
         out = myReduceGWAS(SEL,peakindex);
         out.nPeak = length(peakindex);
         out.nSupport = zeros(1,out.nPeak);
         out.Support = cell(1,out.nPeak);
         for i=1:out.nPeak
            lab = labels(peakindex(i));
            supindex = find(labels==lab);
            if isempty(supindex), continue; end
            out.nSupport(i) = length(supindex);
            out.Support{i}.POS = pos(supindex);
            out.Support{i}.RS = rs(supindex);
         end
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
end



