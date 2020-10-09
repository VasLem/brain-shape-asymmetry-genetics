function out = RemovePeakByMinSupportAndSign(PEAKS,pCrit,minsupport)
         if nargin<3, minsupport=2; end% no LD
         if nargin<2, pCrit = 5e-8;end
         if pCrit<0.5, pCrit = -log10(pCrit);end
         out = [];i = 0;
         for ind=1:1:PEAKS.nPeak
             if PEAKS.nSupport(ind)<minsupport&&PEAKS.bestP(ind)<pCrit, continue; end
             %if PEAKS.bestP(ind)>=pCrit, continue;end
             i = i+1;
             out.RS{i} = PEAKS.RS{ind};
             out.POS(i) = PEAKS.POS(ind);
             out.CHR(i) = PEAKS.CHR(ind);
             out.A1{i} = PEAKS.A1{ind};
             out.A2{i} = PEAKS.A2{ind};
             out.MAF(i) = PEAKS.MAF(ind);
             out.N(i) = PEAKS.N(ind);
             out.bestP(i)= PEAKS.bestP(ind);
             out.P(i,:) = PEAKS.P(ind,:);
             out.SNP(i,:) = PEAKS.SNP(ind,:);
             out.Support{i} = PEAKS.Support{ind};
             out.nSupport(i) = PEAKS.nSupport(ind);
         end
         out.nPeak = i;
end

             