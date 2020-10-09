function out = ExtractPeakByLocationAndPCorrelation(RS,CHR,POS,P,pCrit,pSUG,dist,cT)
         if nargin<8, cT = 0; end% no phenotypic correlation into account
         if nargin<7, dist = 250e3;end
         if nargin<6, pSUG = 5e-7;end
         if nargin<5, pCrit = 5e-8;end
         if pCrit<1, pCrit = -log10(pCrit);end
         if pSUG<1, pSUG = -log10(pSUG);end
         
         bestP = max(P,[],2);
         % First Filtering
         index = find(bestP>=pSUG);
         rs = RS(index);
         chr = CHR(index);
         pos = POS(index);
         p = P(index,:);
         bestp = bestP(index);
         % initiate
         nHits = length(index);
         list = 1:nHits;counter = 0;lab = 0;
         labels = zeros(1,nHits);
         top = zeros(1,nHits);
         while ~isempty(list)&&lab<=nHits
             counter = counter +1;
             %disp(num2str(counter));
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
             p1 = p(bestind,:);p2 = p(subind,:);% select snps correlated phenotypic effect
             pcorr = zeros(1,length(subind));
             for k=1:length(subind)
                pcorr(k) = corr(p1',p2(k,:)','type','spearman'); 
             end
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