function out = overlapPeaksInOtherGWAS(peaks,GWAS,pCrit,pSUG,dist)
         if nargin<5, dist = 250e3;end
         if nargin<4, pSUG = 5e-7;end
         if nargin<3, pCrit = 5e-8;end
         if pCrit<1, pCrit = -log10(pCrit);end
         if pSUG<1, pSUG = -log10(pSUG);end
         n = length(peaks.POS);
         exactOverlap = nan*zeros(n,4);
         proxOverlap = nan*zeros(n,4);
         proxInfo = cell(n,2);
         exactInfo = cell(n,2);
         for i=1:1:n
             % i=10;
            disp(num2str(i));
            index = find(GWAS.CHR==peaks.CHR(i));
            index2 = find(GWAS.POS==peaks.POS(i));
            index3 = find(GWAS.POS<=(peaks.POS(i)+dist)&GWAS.POS>=(peaks.POS(i)-dist));
            % first exact position
            ind=intersect(index,index2);
            if ~isempty(ind)
               [val,segment] = max(GWAS.P(ind,:));
               exactOverlap(i,1) = val;
               exactOverlap(i,2) = segment;
               exactOverlap(i,3) = val>=pCrit;
               exactOverlap(i,4) = val>=pSUG;
               exactInfo{i,1} = GWAS.RS{ind};
               exactInfo{i,2} = GWAS.POS(ind);
            else
               disp('exact SNP not found');
            end
            ind=intersect(index,index3); 
            if ~isempty(ind)
               [val,~] = max(GWAS.P(ind,:),[],1);
               [val,segment] = max(val);
               proxOverlap(i,1) = val;
               proxOverlap(i,2) = segment;
               proxOverlap(i,3) = val>=pCrit;
               proxOverlap(i,4) = val>=pSUG;
               [val,~] = max(GWAS.P(ind,:),[],2);
               [~,indprox] = max(val);
               proxInfo{i,1} = GWAS.RS{ind(indprox(1))};
               proxInfo{i,2} = GWAS.POS(ind(indprox(1)));
            end 
         end
         out.ExactOverlap = exactOverlap;
         out.ProxOverlap = proxOverlap;
         out.ExactCritSug = ((nansum(exactOverlap(:,3:4),1))/n)*100;
         out.ProxCritSug = ((nansum(proxOverlap(:,3:4),1))/n)*100;
         out.ExactCritSugNr = (nansum(exactOverlap(:,3:4),1));
         out.ProxCritSugNr = (nansum(proxOverlap(:,3:4),1));
         out.proxInfo = proxInfo;
         out.exactInfo = exactInfo;
end

% peaks = PEAKS;