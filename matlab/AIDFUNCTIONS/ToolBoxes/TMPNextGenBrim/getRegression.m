function [M,MI,R2] = getRegression(A,B,Booting)
         if nargin < 3, Booting = (1:size(A,2));end
         [A,B] = eliminateNAN(A,B);
         [~,~,~,~,M,~,~,stats] = plsregress(A,B,size(A,2));
         MI = M(1,:);M = M(1+Booting,:);
         if nargout < 3, return; end
         R2 = zeros(size(M));
         n = size(B,1);
         if size(A,2)>1
             for i=1:1:length(Booting)
                 % building the reduced model
                 index = setdiff(1:size(A,2),Booting(i));
                 AR = A(:,index);
                 [~,~,~,~,~,~,~,statsR] = plsregress(AR,B,size(AR,2));
                 E = statsR.Yresiduals;
                 [~,~,~,~,~,~,~,statsR] = plsregress(AR,A(:,Booting(i)));
                 Afor = statsR.Yresiduals;
                 [~,~,~,~,~,~,~,stats] = plsregress(Afor,E,1);
                 P = E-stats.Yresiduals;
                 avg = repmat(mean(E),n,1);
                 SST = sum((E-avg).^2);
                 SSR = sum((P-avg).^2);
                 R2(i,:) = (SSR./SST);
             end
         else
             P = B-stats.Yresiduals;
             avg = repmat(mean(B),n,1);
             SST = sum((B-avg).^2);
             SSR = sum((P-avg).^2);
             R2 = (SSR./SST);
         end
end