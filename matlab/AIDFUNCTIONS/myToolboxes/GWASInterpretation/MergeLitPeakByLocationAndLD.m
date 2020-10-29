function [out,topindex,labels] = MergeLitPeakByLocationAndLD(PEAKS,dist,ldT)
         if nargin<3, ldT=0.2; end% no LD
         if nargin<2, dist = 2e6;end
         R2 = getLDStructure(PEAKS.SNP);
         % Detecting to merge
         nHits = PEAKS.nPeak;
         list = 1:nHits;counter = 0;lab = 0;
         labels = zeros(1,nHits);top = zeros(1,nHits);
         while ~isempty(list)&&counter<=PEAKS.nPeak
             counter = counter+1;
             lab = lab+1;
             %disp(counter);
             %[~,mind]=max(PEAKS.bestP(list));
             mind = 1;
             bestind = list(mind);
             top(bestind) = 1;
             labels(bestind) = lab;
             bestchr = PEAKS.CHR(bestind);
             subindex = find(PEAKS.CHR(list)==bestchr);% select peaks on the same chromosome
             if isempty(subindex), continue; end % there are no more snps on the same chromosome
             subind = list(subindex);
             bestpos = PEAKS.POS(bestind);
             subpos = PEAKS.POS(subind);
             take = subpos<=(bestpos+dist)&subpos>=(bestpos-dist);% select snps in the same proximity
             if sum(take)==0, continue; end % there are no snps close by
             subind = subind(find(take));
             r2 = R2(bestind,subind);
             subind = subind(find(r2>=ldT));
             if isempty(subind), continue; end
             labels(subind) = lab;
             list = setdiff(list,subind);
         end
         % Merging
         topindex = find(top);% run from start to back, then the order along the genome is preserved.
         out.nPeak = length(topindex);
         for i=1:1:out.nPeak
             % i=1;
             ind = topindex(i);
             out.RS{i} = PEAKS.RS{ind};
             out.POS(i) = PEAKS.POS(ind);
             out.CHR(i) = PEAKS.CHR(ind);
             out.SNP(i,:) = PEAKS.SNP(ind,:);
             if isfield(PEAKS,'Study'), out.Study(i) = PEAKS.Study(ind); end
             %out.JOBID(i) = PEAKS.JOBID(ind);
         end
end
             
