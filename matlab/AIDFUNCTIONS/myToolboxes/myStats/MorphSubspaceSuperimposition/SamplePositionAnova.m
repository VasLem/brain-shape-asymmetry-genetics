function out = SamplePositionAnova(X1,X2,t)
         if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting Centroids
           AvgX1 = mean(X1);
           AvgX2 = mean(X2);
         % F-ratio
           F = getFratio(X1,X2);
         % generating Effect output  
           out.AvgX1 = AvgX1;
           out.AvgX2 = AvgX2;
           out.Difference = AvgX1-AvgX2;
           out.Distance = sqrt(sum(out.Difference.^2));
           out.F = F;
         % Permutation test  
           if t<=0, return; end
           FCount = false(1,t);
           nT = nX1+nX2;
           X = [X1; X2];
           disp('Permuting');
           tic;
           parfor i=1:t
                  ind = randperm(nT);
                  X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
                  X2perm = X(ind(nX1+1:end),:);
                  Fperm = getFratio(X1perm,X2perm);
                  FCount(i) = Fperm>=F;
           end
           toc;
           out.pF = sum(FCount)/t;
end

function D = getDistanceMatrix(X1,X2,type)
         D = squareform(pdist([X1;X2],type));
end


function F = getFratio(X1,X2)
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         N = nX1+nX2;
         D = squareform(pdist([X1;X2],'euclidean'));
         % getting the total SS Following Anderson
         SST = 0;
         for i=1:1:(N-1)
            for j=(i+1):1:N
                SST = SST+D(i,j)^2;
            end
         end
         SST = SST/N;
         % getting the within SS
         tic;
         SSW = 0;
         D = squareform(pdist(X1,'euclidean'));
         for i=1:1:(nX1-1)
            for j=(i+1):1:nX1
                SSW = SSW+D(i,j)^2;
            end
         end 
         D = squareform(pdist(X2,'euclidean'));
         for i=1:1:(nX2-1)
            for j=(i+1):1:nX2
                SSW = SSW+D(i,j)^2;
            end
         end
         SSW = SSW/nX1;
         toc;
%          SSW = 0;
%          for i=1:1:(N-1)
%             for j=(i+1):1:N
%                 if i<=nX1&&j<=nX1
%                    e = 1;
%                 elseif i>=nX1&&j>=nX1
%                    e = 1;
%                 else
%                    e = 0;
%                 end
%                 SSW = SSW+e*D(i,j)^2;
%             end
%          end
%          SSW = SSW/nX1;
         % Getting F-ratio
         F = (SST-SSW)/(SSW/(N-2));
end

% function out = SubspacePosition(X1,X2,type,t)
%          if nargin < 4, t = 0; end
%          nX1 = size(X1,1);
%          nX2 = size(X2,1);
%          nV = size(X2,2);
%          % Getting Subspace Resroids & Residus
%              AvgX1 = mean(X1);
%              AvgX2 = mean(X2);
%              ResX1 = X1-repmat(AvgX1,nX1,1);
%              ResX2 = X2-repmat(AvgX2,nX2,1);
%          % Getting within Subspace Covariances + pooled 
%              switch lower(type)
%                  case 'euclidean'
%                      S = eye(nV,nV);
%                  case 'normalized euclidean'
%                      W1 = diag(std(ResX1));
%                      W2 = diag(std(ResX2));
%                      S = inv((nX1*W1 + nX2*W2)/(nX1+nX2));
%                  case 'mahalanobis'
%                      W1 = cov(ResX1);
%                      W2 = cov(ResX2);
%                      S = inv((nX1*W1 + nX2*W2)/(nX1+nX2));
%              end
%          % Getting distance between averages
%            EffectSize = sqrt((AvgX1-AvgX2)*S*(AvgX1-AvgX2)');
%          % generating Effect output  
%            out.AvgX1 = AvgX1;
%            out.AvgX2 = AvgX2;
%            out.Difference = AvgX1-AvgX2;
%            out.EffectSize = EffectSize;
%          % Permutation test  
%            if t<=0, return; end
%            % generating test-statistic (Euclidean Distance)
%            DTest = sqrt((AvgX1-AvgX2)*(AvgX1-AvgX2)'); 
%            TestCount = false(1,t);
%            nT = nX1+nX2;
%            X = [X1; X2];
%            disp('Permuting');
%            tic;
%            parfor i=1:t
%                   ind = randperm(nT);
%                   X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
%                   X2perm = X(ind(nX1+1:end),:);
%                   AvgX1perm = mean(X1perm);
%                   AvgX2perm = mean(X2perm);
%                   Dperm = sqrt((AvgX1perm-AvgX2perm)*(AvgX1perm-AvgX2perm)');
%                   TestCount(i) = Dperm>=DTest;
%            end
%            toc;
%            out.DTest = DTest;
%            out.pD = sum(TestCount)/t;
% end