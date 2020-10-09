function [stat,se,p] = myLDSCOV(Z1,Z2,L2,n)
       
       Z1 = Z1(:);
       Z2 = Z2(:);
       L2= L2(:);

       nSNP = length(L2);
       %W = 1./L2(:);
       W = ones(nSNP,1);
       X = lscov([ones(nSNP,1) L2 Z1],Z2,W);
       X1 = X(end);
       X = lscov([ones(nSNP,1) L2 Z2],Z1,W);
       X2 = X(end);
       stat = (X1+X2)/2; 
       %n = 1000;
       val = zeros(1,n);
       parfor i=1:n
           index = randsample(nSNP,nSNP,true);
           X = lscov([ones(nSNP,1) L2(index) Z1(index)],Z2(index),W(index));
           X1 = X(end);
           X = lscov([ones(nSNP,1) L2(index) Z2(index)],Z1(index),W(index));
           X2 = X(end);
           val(i) = (X1+X2)/2;   
       end
       se = std(val);
       p = 2*normcdf(-1*abs(stat/se),0,1);
end


% function [stat,se,p] = myLDSCOV(Z1,Z2,L2,LDblockID,Blocks)
%        Z1 = Z1(:);
%        Z2 = Z2(:);
%        L2= L2(:);
%        nSNP = length(L2);
%        W = ones(nSNP,1);
%        X = lscov([ones(nSNP,1) L2 Z1],Z2,W);
%        X1 = X(end);
%        X = lscov([ones(nSNP,1) L2 Z2],Z1,W);
%        X2 = X(end);
%        stat = (X1+X2)/2; 
%        n = length(Blocks);
%        val = zeros(1,n);
%        parfor i=1:n
%            index = find(~(LDblockID==Blocks(i)));
%            X = lscov([ones(length(index),1) L2(index) Z1(index)],Z2(index),W(index));
%            X1 = X(end);
%            X = lscov([ones(length(index),1) L2(index) Z2(index)],Z1(index),W(index));
%            X2 = X(end);
%            val(i) = (X1+X2)/2;   
%        end
%        se = std(val);
%        p = 2*normcdf(-1*abs(stat/se),0,1);
% end
