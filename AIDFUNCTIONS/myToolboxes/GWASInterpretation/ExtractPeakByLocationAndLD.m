function out = ExtractPeakByLocationAndLDAndPCorr(SEL,pCrit,dist,ldT)
         if nargin<4, ldT=0; end% no LD
         if nargin<3, dist = 500e3;end
         if nargin<2, pCrit = 5e-8;end
         if pCrit<1, pCrit = -log10(pCrit);end
      
         rs = SEL.RS;
         chr = SEL.CHR;
         pos = SEL.POS;
         p = SEL.P;
         bestp = SEL.bestP;
         SNP = SEL.SNP;
         SNP_1KG = SEL.SNP1KG;
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
             take = subpos<=(bestpos+dist)&subpos>=(bestpos-dist);% select snps in the same proximity
             if sum(take)==0, continue; end % there are no snps close by
             subind = subind(find(take));
             
             
             GENO1 = SNP(bestind,:);
             GENO2 = SNP(subind,:);
             R2 = getLDStructure(GENO2,GENO1);
             
             GENO1_1KG = SNP_1KG(bestind,:);
             GENO2_1KG = SNP_1KG(subind,:);
             R2_1KG = getLDStructure(GENO2_1KG,GENO1_1KG);
             
             
             GENO1 = GENO2_1KG;
             GENO2 = GENO1_1KG;
             
             GENO1 = GENO2_1KG;
             GENO2 = GENO2_1KG;
             
             
             R2TOT = getLDStructure(GENO2);
             
             R2TOT_1KG = getLDStructure(GENO2_1KG);
             
             
             figure;plot(1:length(subind),R2,'b.','MarkerSize',15);
             
             
             figure;scatter(1:length(subind),sqrt(R2),30,bestp(subind),'filled');set(gca,'clim',[7.5 30]);colormap('jet');
             
             
             figure;scatter(1:length(subind),bestp(subind),30,R2,'filled');set(gca,'ylim',[5 30]);colormap('jet');
             
             figure;scatter(1:length(subind),bestp(subind),30,pcorr,'filled');set(gca,'ylim',[5 30]);colormap('jet');
             
             
             
             p1 = p(bestind,:);p2 = p(subind,:);% select snps correlated phenotypic effect
             pcorr = zeros(1,length(subind));
             for k=1:length(subind)
                pcorr(k) = corr(p1',p2(k,:)','type','spearman'); 
             end
             
             
             pcorr = zeros(length(subind),length(subind));
             parfor k1=1:length(subind)
                 tmp = zeros(1,length(subind));
                 for k2=1:length(subind)
                    tmp(k2) = corr(p2(k1,:)',p2(k2,:)','type','spearman'); 
                 end
                 pcorr(k1,:) = tmp;
             end
             
             
             figure;plot(pcorr(:),sqrt(R2TOT(:)),'b.','MarkerSize',15);
             
             
             
             subind = subind(find(pcorr>=cT));
             labels(subind) = lab;
             list = setdiff(list,subind);
         end
         peakindex = find(top);
         subindex = find(bestp(peakindex,:)>=pCrit);
         peakindex = peakindex(subindex); % only retain peaks with pCRIT top
         out.nPeak = length(peakindex);
         out.RS = rs(peakindex);
         out.POS = pos(peakindex);
         out.CHR = chr(peakindex);
         out.bestP = bestp(peakindex);
         out.P = p(peakindex,:);
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