function [out] = Match2GIMP(input,fac)
         if nargin<2, fac = 5; end
         nF = size(input,2);out.nF = nF;
         nS = size(input,1);out.nS = nS;
         
         ID = eye(nF);
         indgen = find(ID==1);ngen = length(indgen);out.nGEN = ngen;
         indimpos = find(ID==0);nimpos = length(indimpos);out.nIMPOS = nimpos;
         GEN = zeros(nS,ngen);
         IMPOS = zeros(nS,nimpos);
         for s=1:1:nS
             tmp = squeeze(input(s,:,:));
             GEN(s,:) = tmp(indgen);
             IMPOS(s,:) = tmp(indimpos);
         end
         out.GEN = GEN;
         out.IMPOS = IMPOS;
         out.REDIMPOS = IMPOS;
%          out.REDIMPOS = setdiff(IMPOS',GEN','rows')';
         if size(out.REDIMPOS,2)>fac*ngen
            rng(3);ind = randsample(size(out.REDIMPOS,2),ngen*fac); 
            out.REDIMPOS = out.REDIMPOS(:,ind); 
         end
end