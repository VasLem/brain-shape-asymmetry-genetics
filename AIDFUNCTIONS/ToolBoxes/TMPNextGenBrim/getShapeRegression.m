function [M,R2] = getShapeRegression(A,B,Booting)
         if nargin < 3, Booting = (1:size(A,2));end
         [A,B] = eliminateNAN(A,B);
         [~,~,~,~,M] = plsregress(A,B,size(A,2));
         M = M(1+Booting,:);
         if nargout < 2, return; end
         R2 = zeros(size(M,1),size(M,2)/3);
         for i=1:1:length(Booting)
             % building the reduced model
             index = setdiff(1:size(A,2),Booting(i));
             AR = A(:,index);
             [~,~,~,~,~,~,~,statsR] = plsregress(AR,B,size(AR,2));
             E = statsR.Yresiduals;
             [~,~,~,~,~,~,~,statsR] = plsregress(AR,A(:,Booting(i)));
             Afor = statsR.Yresiduals;
             [~,~,~,~,~,~,~,stats] = plsregress(Afor,E,1);
             [n,nE] = size(E);
             P = E-stats.Yresiduals;
             avg = repmat(mean(E),n,1);
             SST = sum((E-avg).^2);
             SSR = sum((P-avg).^2);
             SSR = reshape(SSR,3,nE/3);
             SSR = sum(SSR);
             SST = reshape(SST,3,nE/3);
             SST = sum(SST);
             R2(i,:) = (SSR./SST);
         end
end