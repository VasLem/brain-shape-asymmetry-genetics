function out = ExtractPeakByLocation(RS,CHR,POS,P,pCrit,dist)
         if nargin<6, dist = 250e3;end
         if nargin<5, pCrit = 5e-8;end
         if pCrit<1, pCrit = -log10(pCrit);end
         if size(P,2)>1, P = max(P,[],2); end % Given pvalues are already -log transformed
         % First Filtering
         index = find(P>=pCrit);
         rs = RS(index);
         chr = CHR(index);
         pos = POS(index);
         p = P(index);
         
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
             [m,mind] = max(p(list));
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
             labels(subind) = lab;
             list = setdiff(list,subind);
         end
         out.nPeak = sum(top);
         peakindex = find(top);
         out.RS = rs(peakindex);
         out.POS = pos(peakindex);
         out.CHR = chr(peakindex);
         out.P = p(peakindex);     
end