function D = GTDistance(v1,v2,fAA,faA,faa)      
         if nargin < 3
             freq = ones(size(v1));
         else
             freq = nan*zeros(size(v1));
             freq(v1==-1) = 1-fAA(v1==-1);
             freq(v1==0) = 1-faA(v1==0);
             freq(v1==1) = 1-faa(v1==1);
         end
         s = size(v2,1);
         ind = find(~isnan(v1));
         v1 = v1(ind);
         v2 = v2(:,ind);
         freq = freq(ind);
         N = repmat(freq,s,1);
         N(isnan(v2)) = nan;
         tmp = (repmat(v1,s,1)==v2).*repmat(freq,s,1);
         tmp(isnan(v2)) = nan;
         D = 1-nansum(tmp,2)./nansum(N,2);
end