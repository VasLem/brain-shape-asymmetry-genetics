function [GC,pGC,seGC] = geneticCorrelationBootstrapLDSCOV(V1,V2,L2,n)
         nV1 = size(V1,2);
         nV2 = size(V2,2);
         GC = zeros(nV1,nV2);
         pGC = ones(nV1,nV2);
         seGC = ones(nV1,nV2);
         disp('COMPUTING PAIRWISE CORRELATIONS');
         [path,ID] = setupParForProgress(nV1);
         parfor i=1:nV1
             % i=1;
             vv1 = V1(:,i);
             forGC = zeros(1,nV2);
             forpGC = ones(1,nV2);
             forseGC = ones(1,nV2);
             for j=1:nV2
                 % j=7; 
                 vv2 = V2(:,j); %#ok<*PFBNS>
                 [forGC(j),forseGC(j),forpGC(j)] = innerLDSCOV(vv1(:),vv2(:),L2,n);
             end
             GC(i,:) = forGC;
             pGC(i,:) = forpGC;
             seGC(i,:) = forseGC;
             parfor_progress;
         end
         closeParForProgress(path,ID)
end

function [stat,se,p] = innerLDSCOV(vv1,vv2,L2,K)   
       %vv1 = myGetRank(vv1(:));
       %vv2 = myGetRank(vv2(:));     
       vv1 = vv1(:);
       vv2 = vv2(:);
       
       L2= L2(:);
       X = lscov([ones(length(L2),1) L2 vv1],vv2);
       X1 = X(end);
       X = lscov([ones(length(L2),1) L2 vv2],vv1);
       X2 = X(end);
       stat = (X1+X2)/2;
       n = length(vv1);
       nPerK = round(n/K);
       ind = 1:nPerK:n;
       val = zeros(1,K);
       parfor i=1:K-1
           %i=1
          index = setdiff(1:n,ind(i):ind(i+1)); 
          tic;X = lscov([ones(length(index),1) L2(index) vv1(index)],vv2(index));
          X1 = X(end);
          X = lscov([ones(length(index),1) L2(index) vv2(index)],vv1(index));
          X2 = X(end);
          val(i) = (X1+X2)/2;
       end
       se = std(val);
       p = 2*normcdf(-1*abs(stat/se),0,1);
end

function r = myGetRank(Data)
         data_sorted = sort(Data);
         [~, rnk] = ismember(Data,data_sorted);
         r = rnk(:);
end


%        n = length(vv1);
%        K = 10;
%        c = cvpartition(n,'KFold',K); 
%        val = zeros(1,K);
%        parfor i=1:K
%            index = find(training(c,i));
%            nSNP = length(index);
%            X = lscov([ones(nSNP,1) L2(index) vv1(index)],vv2(index),ones(nSNP,1));
%            X1 = X(end);
%            X = lscov([ones(nSNP,1) L2(index) vv2(index)],vv1(index),ones(nSNP,1));
%            X2 = X(end);
%            val(i) = (X1+X2)/2; 
%        end
%        se = std(val);
%        %p = 2*normcdf(-1*abs(stat/se),0,1);
%        p = normcdf(-1*(stat/se),0,1);


%        n = 100;
%        val = zeros(1,n);
%        parfor i=1:n
%            index = randsample(nSNP,nSNP,true);
%            X = lscov([ones(nSNP,1) L2(index) vv1(index)],vv2(index),W(index));
%            X1 = X(end);
%            X = lscov([ones(nSNP,1) L2(index) vv2(index)],vv1(index),W(index));
%            X2 = X(end);
%            val(i) = (X1+X2)/2;   
%        end
%        se = std(val);
%        p = 2*normcdf(-1*abs(stat/se),0,1);
% %       se = 0;
% %       p = 0;


% function [stat,se,p] = innerLDSCOV(Z1,Z2,L2,n)   
%        Z1 = Z1(:);
%        Z2 = Z2(:);
%        L2= L2(:);
%        nSNP = length(L2);
%        %W = 1./L2(:);
%        W = ones(nSNP,1);
%        X = lscov([ones(nSNP,1) L2 Z1],Z2,W);
%        X1 = X(end);
%        X = lscov([ones(nSNP,1) L2 Z2],Z1,W);
%        X2 = X(end);
%        stat = (X1+X2)/2; 
%        %n = 1000;
% %        val = zeros(1,n);
% % %        parfor i=1:n
% % %            index = randsample(nSNP,nSNP,true);
% % %            X = lscov([ones(nSNP,1) L2(index) Z1(index)],Z2(index),W(index));
% % %            X1 = X(end);
% % %            X = lscov([ones(nSNP,1) L2(index) Z2(index)],Z1(index),W(index));
% % %            X2 = X(end);
% % %            val(i) = (X1+X2)/2;   
% % %        end
% %        se = std(val);
% %        p = 2*normcdf(-1*abs(stat/se),0,1);
%       se = 0;
%       p = 0;
% end

% function [stat,se,p] = getJacknifeSTATS(vv1,vv2)
%          stat = corr(vv1(:),vv2(:),'type','Spearman');
%          n = length(vv1);
%          K = 20;
%          c = cvpartition(n,'KFold',K); 
%          val = zeros(1,K);
%          for i=1:K
%              index = training(c,i);
%              val(i)= corr(vv1(index),vv2(index),'type','Spearman');
%          end
%          se = std(val);
%          %p = 2*normcdf(-1*abs(stat/se),0,1);
%          p = normcdf(-1*(stat/se),0,1);
% end

% function [stat,se,p] = myLDSCOVwithSE(Z1,Z2,L2,LDblockID,Blocks)
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