function [TMatches,FMatches] = seperateMatches(Matches)
         [nM,n,~] = size(Matches);
         ind1 = repmat((1:nM)',1,n);
         ind2 = repmat((1:n),nM,1);
         Aind = 1:(n*n*nM);
         Tind = sub2ind([nM n n],ind1(:),ind2(:),ind2(:));
         Find = setdiff(Aind,Tind);
         TMatches = reshape(Matches(Tind),nM,n);
         FMatches = reshape(Matches(Find),nM,n,n-1);
end