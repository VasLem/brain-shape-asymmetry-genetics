function out = SubspacePositionDistance(X1,X2,type,t)
         if nargin < 4, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         nV = size(X2,2);
         % Getting Subspace Resroids & Residus
             AvgX1 = mean(X1);
             AvgX2 = mean(X2);
             ResX1 = X1-repmat(AvgX1,nX1,1);
             ResX2 = X2-repmat(AvgX2,nX2,1);
         % Getting within Subspace Covariances + pooled 
             switch lower(type)
                 case 'euclidean'
                     S = eye(nV,nV);
                 case 'normalized euclidean'
                     W1 = diag(std(ResX1));
                     W2 = diag(std(ResX2));
                     S = inv((nX1*W1 + nX2*W2)/(nX1+nX2));
                 case 'mahalanobis'
                     W1 = cov(ResX1);
                     W2 = cov(ResX2);
                     S = inv((nX1*W1 + nX2*W2)/(nX1+nX2));
             end
         % Getting distance between averages
           EffectSize = sqrt((AvgX1-AvgX2)*S*(AvgX1-AvgX2)');
         % generating Effect output  
           out.AvgX1 = AvgX1;
           out.AvgX2 = AvgX2;
           out.Difference = AvgX1-AvgX2;
           out.EffectSize = EffectSize;
         % Permutation test  
           if t<=0, return; end
           % generating test-statistic (Euclidean Distance)
           DTest = sqrt((AvgX1-AvgX2)*(AvgX1-AvgX2)'); 
           TestCount = false(1,t);
           nT = nX1+nX2;
           X = [X1; X2];
           disp('Permuting');
           tic;
           parfor i=1:t
                  ind = randperm(nT);
                  X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
                  X2perm = X(ind(nX1+1:end),:);
                  AvgX1perm = mean(X1perm);
                  AvgX2perm = mean(X2perm);
                  Dperm = sqrt((AvgX1perm-AvgX2perm)*(AvgX1perm-AvgX2perm)');
                  TestCount(i) = Dperm>=DTest;
           end
           toc;
           out.DTest = DTest;
           out.pD = sum(TestCount)/t;
end