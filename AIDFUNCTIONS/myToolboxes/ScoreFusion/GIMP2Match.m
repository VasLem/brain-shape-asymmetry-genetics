function out = GIMP2Match(input)
         nF = input.nF;
         nS = input.nS;
         ID = eye(nF);
         indgen = find(ID==1);ngen = length(indgen);out.nGEN = ngen;
         indimpos = find(ID==0);nimpos = length(indimpos);out.nIMPOS = nimpos;
         
         out = nan*zeros(nS,nF,nF);
         for s=1:1:nS
             tmp = nan*zeros(nF,nF);
             tmp(indgen) = input.GEN(s,:);
             tmp(indimpos) = input.IMPOS(s,:);
             out(s,:,:) = tmp;
         end
end