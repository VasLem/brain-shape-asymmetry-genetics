% Katleen: MixedMomentDStatistic van Peter veranderd zodat hij Paired
% werkt
function out = MixedMomentDStatisticPaired(X1,X2,t,mcp,corr)
         if nargin < 5, corr = false; end
         if nargin < 4, mcp = +inf; end
         if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         mcp = min([nX1 nX2 mcp]);
         % Standarizing the Data
           % Centering
%            disp('centering');
           X1 = X1-repmat(mean(X1),nX1,1);
           X2 = X2-repmat(mean(X2),nX2,1);
           if corr% scaling
              disp('scaling');
              X1 = X1./repmat(std(X1),nX1,1);
              X2 = X2./repmat(std(X2),nX2,1);
           end
         % Create Subspaces
           [SSpaceX1,EigValX1] = createSubspaceSVD(X1);
           [SSpaceX2,EigValX2] = createSubspaceSVD(X2);
           SSpaceX1 = SSpaceX1(:,1:mcp);EigValX1 = EigValX1(1:mcp);
           SSpaceX2 = SSpaceX2(:,1:mcp);EigValX2 = EigValX2(1:mcp);
         % Getting Distance Statistic
           [DStat,A] = getSubspaceDistances(SSpaceX1,SSpaceX2,'projection');
           EigValDStat = cumsum(abs(EigValX1-EigValX2));
         % generating Effect output  
           out.SSpaceX1 = SSpaceX1;
           out.EigValX1 = EigValX1;
           out.SSpaceX2 = SSpaceX2;
           out.EigValX2 = EigValX2;
           out.DStat = DStat;
           out.A = A;
           out.EigValDStat = EigValDStat;
         % Permutation test  
           if t<=0, return; end
           StatCount = false(length(DStat),t);
           EigValStatCount = false(length(EigValDStat),t);
           DStatpermC = zeros(t,length(DStat));
           EigValDStatpermC = zeros(t,length(EigValDStat));
           nT = nX1+nX2;
           X = [X1; X2];
           disp('Permuting');
           tic;
           parfor i=1:t
                  r = randi(2,nX1,1);
                  X1perm = zeros(size(X1));
                  X2perm = zeros(size(X2));
                  for k=1:1:nX1
                      switch r(k)
                          case 1
                              X1perm(k,:) = X1(k,:);
                              X2perm(k,:) = X2(k,:);
                          case 2
                              X1perm(k,:) = X2(k,:);
                              X2perm(k,:) = X1(k,:);
                      end
                  end
                  [SSpaceX1perm,EigValX1perm] = createSubspaceSVD(X1perm);
                  [SSpaceX2perm,EigValX2perm] = createSubspaceSVD(X2perm);
                  SSpaceX1perm = SSpaceX1perm(:,1:mcp);EigValX1perm = EigValX1perm(1:mcp);
                  SSpaceX2perm = SSpaceX2perm(:,1:mcp);EigValX2perm = EigValX2perm(1:mcp);
                  DStatperm = getSubspaceDistances(SSpaceX1perm,SSpaceX2perm,'projection');
                  EigValDStatperm = cumsum(abs(EigValX1perm-EigValX2perm));
                  StatCount(:,i) = (DStatperm>=DStat)';
                  EigValStatCount(:,i) = EigValDStatperm>=EigValDStat;
                  DStatpermC(i,:) = DStatperm;
                  EigValDStatpermC(i,:) = EigValDStatperm;
           end
           toc;
           figure;
           %subplot(1,2,1);
           set(gca,'Xlim',[1 mcp]);
           hold on;
           for i=1:1:t
               plot(DStatpermC(i,:),'r-.');
               drawnow;   
           end
           plot(DStat,'b','linewidth',3);
           %subplot(1,2,2);
           %set(gca,'Xlim',[1 mcp]);
           %hold on;
           %for i=1:1:t
           %    plot(EigValDStatpermC(i,:),'r-.');
           %    drawnow;   
           %end
           %plot(EigValDStat,'b','linewidth',3);
           out.pperm = (sum(StatCount,2)/t)';
           out.EigValpperm = (sum(EigValStatCount,2)/t)';
           out.DStatpermC = DStatpermC;
           out.EigValDStatpermC = EigValDStatpermC;
end