function out = SampleSizeAnova(X1,X2,t)
        if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting Subspace centroids & residus
             AvgX1 = mean(X1);
             AvgX2 = mean(X2);
             ResX1 = X1-repmat(AvgX1,nX1,1);
             ResX2 = X2-repmat(AvgX2,nX2,1);
         % Getting within group distances
             PDX1 = sqrt(sum(ResX1.^2,2));
             PDX2 = sqrt(sum(ResX2.^2,2));
         % Getting mean Distances
             out.DispX1 = mean(PDX1);
             out.DispX2 = mean(PDX2);
             X = [PDX1(:);PDX2(:)];
             G = [ones(nX1,1);zeros(nX2,1)];
             [~,TAB] = anova1(X,G,'off');
             F = TAB{2,5};
             out.F = F;
             if t<=0, return; end
         % permuting
             FCount = false(1,t);
             nT = nX1+nX2;
             disp('Permuting');
             tic;
             parfor i=1:t
                 ind = randperm(nT);
                 [~,TABfor] = anova1(X(ind),G,'off'); %#ok<*PFBNS>
                 FCount(i) = TABfor{2,5} >= F;
             end
             toc;
             out.pF = (sum(FCount)/t);
end
% function out = SubspaceSize(X1,X2,t)
%          type = 'euclidean';
%          if nargin < 3, t = 0; end
%          nX1 = size(X1,1);
%          nX2 = size(X2,1);
%          nV = size(X2,2);
%          % Getting Subspace Resroids
%              AvgX1 = mean(X1);
%              AvgX2 = mean(X2);
%              ResX1 = X1-repmat(AvgX1,nX1,1);
%              ResX2 = X2-repmat(AvgX2,nX2,1);
%          % Getting within Subspace Covariances + pooled 
%              switch lower(type)
%                  case 'euclidean'
%                      S = eye(nV,nV);
%                      W1 = eye(nV,nV);
%                      W2 = eye(nV,nV);
%                  case 'normalized euclidean'% This makes no sense to do
%                      W1 = inv(diag(std(ResX1)));
%                      W2 = inv(diag(std(ResX2)));
%                  case 'mahalanobis'% This makes no sense to do
%                      W1 = inv(cov(ResX1));
%                      W2 = inv(cov(ResX2));
%              end
%          % Getting within group distances
%              PDX1 = zeros(1,nX1);
%              parfor i=1:nX1
%                  PDX1(i) = sqrt(ResX1(i,:)*W1*ResX1(i,:)');
%              end
%              PDX2 = zeros(1,nX2);
%              parfor i=1:nX2
%                  PDX2(i) = sqrt(ResX2(i,:)*W2*ResX2(i,:)');
%              end
%          % Getting mean Distances
%              out.DispX1 = mean(PDX1);
%              out.DispX2 = mean(PDX2);
%              X = [PDX1(:);PDX2(:)];
%              G = [ones(nX1,1);zeros(nX2,1)];
%              [~,TAB] = anova1(X,G,'off');
%              F = TAB{2,5};
%              out.F = F;
%              if t<=0, return; end
%          % permuting
%              FCount = false(1,t);
%              nT = nX1+nX2;
%              disp('Permuting');
%              tic;
%              parfor i=1:t
%                  ind = randperm(nT);
%                  [~,TABfor] = anova1(X(ind),G,'off'); %#ok<*PFBNS>
%                  FCount(i) = TABfor{2,5} >= F;
%              end
%              toc;
%              out.pF = (sum(FCount)/t);
% end
%                           