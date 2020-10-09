function out = SubspaceOrientation(X1,X2,type,t,ncomp)
         if nargin < 5, ncomp = +inf; end
         if nargin < 4, t = 0; end
         %type = 'Procrustes';
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting Subspace centroids & Residus
             AvgX1 = mean(X1);
             AvgX2 = mean(X2);
             ResX1 = X1-repmat(AvgX1,nX1,1);
             ResX2 = X2-repmat(AvgX2,nX2,1);
         % Getting Subspace Scale
             SX1 = mean(sqrt(sum(ResX1.^2,2)));
             SX2 = mean(sqrt(sum(ResX2.^2,2)));
             ResX1 = ResX1/SX1;
             ResX2 = ResX2/SX2;
         % Getting principal angles
             A = my_subspacea(ResX1',ResX2',ncomp);
             D = getDistance(A,type);
             out.A = A;
             out.D = D;
         % permuting
             if t<=0, return; end
             DCount = false(1,t);
             nT = nX1+nX2;
             X = [ResX1; ResX2];
             disp('Permuting');
             tic;
             parfor i=1:t
                  ind = randperm(nT);
                  ResX1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
                  ResX2perm = X(ind(nX1+1:end),:);
                  Aperm = my_subspacea(ResX1perm',ResX2perm',ncomp);
                  Dperm = getDistance(Aperm,type);
                  DCount(i) = Dperm>=D;
             end
             toc;
             out.pD = sum(DCount)/t;
end


function D = getDistance(A,type)
     switch lower(type)
         case 'projection'
             %D = ProjectionMetric(A);
             D = sqrt(length(A)-sum(cos(A).^2));
         case 'procrustes'
             %D = ProcrustesMetric(A);
             D = 2*sqrt(sum(sin(A/2).^2));
         case 'mitteroecker'
             %D = MitteroeckerMetric(A);
             D = sqrt(sum(log(cos(A)).^2));
         case 'average'
             D = length(A)/sum(cos(A));
     end
end

% function out = SubspaceOrientation(X1,X2,type,t,ncomp)
%          if nargin < 5, ncomp = +inf; end
%          if nargin < 4, t = 0; end
%          %type = 'Procrustes';
%          nX1 = size(X1,1);
%          nX2 = size(X2,1);
%          nV = size(X1,2);
%          % Getting Subspace centroids & Residus
%              AvgX1 = mean(X1);
%              AvgX2 = mean(X2);
%              ResX1 = X1-repmat(AvgX1,nX1,1);
%              ResX2 = X2-repmat(AvgX2,nX2,1);
%          % Getting Subspace Scale
%              SX1 = mean(sqrt(sum(ResX1.^2,2)));
%              SX2 = mean(sqrt(sum(ResX2.^2,2)));
%              ResX1 = ResX1/SX1;
%              ResX2 = ResX2/SX2;
%          % Getting principal angles
%              A = my_subspacea(ResX1',ResX2',ncomp);
%              D = getDistance(A,type);
%              out.A = A;
%              out.D = D;
%          % Getting Random Angles
%              dim = nV;
%              tR = 10000;
%              RandAngles = zeros(1,tR);
%              parfor i=1:tR
%                    Dir1 = -1 + 2*rand(1,dim);Dir1 = Dir1/norm(Dir1);
%                    Dir2 = -1 + 2*rand(1,dim);Dir2 = Dir2/norm(Dir2);
%                    RandAngles(i) = angle(Dir1',Dir2');
%              end
%              Pvalues = zeros(1,length(A));
%              parfor i=1:length(Pvalues)
%                     Count = abs(RandAngles)>=cos(A(i));
%                     Pvalues(i) = sum(Count)/tR;
%              end
%              Rank = sum(Pvalues<=0.001);
%          % permuting
%              if t<=0, return; end
%              DCount = false(1,t);
%              nT = nX1+nX2;
%              X = [ResX1; ResX2];
%              disp('Permuting');
%              tic;
%              parfor i=1:t
%                   ind = randperm(nT);
%                   ResX1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
%                   ResX2perm = X(ind(nX1+1:end),:);
%                   Aperm = my_subspacea(ResX1perm',ResX2perm',ncomp);
%                   Dperm = getDistance(Aperm,type);
%                   DCount(i) = Dperm>=D;
%              end
%              toc;
%              out.pD = sum(DCount)/t;
% end
% 
% 
% function D = getDistance(A,type)
%      switch lower(type)
%          case 'projection'
%              D = ProjectionMetric(A);
%          case 'procrustes'
%              D = ProcrustesMetric(A);
%          case 'mitteroecker'
%              D = MitteroeckerMetric(A);
%      end
% end

