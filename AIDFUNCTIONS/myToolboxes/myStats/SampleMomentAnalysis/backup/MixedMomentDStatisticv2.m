function out = MixedMomentDStatisticv2(X1,X2,t,mcp)
         if nargin < 4, mcp = +inf; end
         if nargin < 3, t = 0; end
         type = 'procrustes';
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         mcp = min([nX1 nX2 mcp]);
         % Create Subspaces
           SSpaceX1 = createSubspace(X1');
           SSpaceX2 = createSubspace(X2');
           SSpaceX1 = SSpaceX1(:,1:mcp);
           SSpaceX2 = SSpaceX2(:,1:mcp);
         % Getting Distance Statistic
           [DStat,A] = getSubspaceDistances(SSpaceX1,SSpaceX2,type);
         % generating Effect output  
           out.SSpaceX1 = SSpaceX1;
           out.SSpaceX2 = SSpaceX2;
           out.Dstat = DStat;
           out.A = A;
         % Permutation test  
           if t<=0, return; end
           StatCount = false(length(DStat),t);
           DStatpermC = zeros(t,length(DStat));
           nT = nX1+nX2;
           X = [X1; X2];
           disp('Permuting');
           tic;
           parfor i=1:t
                  ind = randperm(nT);
                  X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
                  X2perm = X(ind(nX1+1:end),:);
                  SSpaceX1perm = createSubspace(X1perm');
                  SSpaceX2perm = createSubspace(X2perm');
                  SSpaceX1 = SSpaceX1(:,1:mcp);
                  SSpaceX2 = SSpaceX2(:,1:mcp);
                  DStatperm = getSubspaceDistances(SSpaceX1perm,SSpaceX2perm,type);
                  StatCounttmp = DStatperm>=DStat;
                  StatCount(:,i) = StatCounttmp';
                  DStatpermC(i,:) = DStatperm;
           end
           toc;
           figure;plot(DStat,'b','linewidth',3);
           hold on;
           for i=1:1:t
               plot(DStatpermC(i,:),'r-.');
               drawnow;   
           end
           out.pperm = (sum(StatCount,2)/t)';
end