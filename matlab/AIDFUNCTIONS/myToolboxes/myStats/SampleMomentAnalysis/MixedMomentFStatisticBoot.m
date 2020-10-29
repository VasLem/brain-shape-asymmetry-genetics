function out = MixedMomentFStatisticBoot(X1,X2,tb,t,ncmp,corr)
         if nargin < 4, t = 0; end
         if nargin < 5, ncmp = +inf; end
         if nargin < 6, corr = false; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         ncmp = min([nX1 nX2 ncmp]);
         % Standarizing the Data
           % Centering
           disp('centering');
           X1 = X1-repmat(mean(X1),nX1,1);
           X2 = X2-repmat(mean(X2),nX2,1);
           if corr% scaling
              disp('scaling');
              X1 = X1./repmat(std(X1),nX1,1);
              X2 = X2./repmat(std(X2),nX2,1);
           end
         % Bootstrapping Subspaces
           SpacesX1 = cell(1,tb);
           SpacesX2 = cell(1,tb);
           disp('Booting Subspaces');
           tic;
           parfor i=1:tb
               ind = randsample((1:nX1),nX1,true);
               tmp = createSubspaceSVD(X1(ind,:)); %#ok<*PFBNS>
               SpacesX1{i} = tmp(:,1:ncmp);
               ind = randsample((1:nX2),nX2,true);
               tmp = createSubspaceSVD(X2(ind,:));       
               SpacesX2{i} = tmp(:,1:ncmp);
           end
           toc;
         % Generating the DistanceMatrix
             Spaces = [SpacesX1';SpacesX2'];
             N = size(Spaces,1);
             DU = zeros(N,N,ncmp);
             disp('Computing Distance Matrix');
             tic;
             parfor i=1:N-1
                 Drow = zeros(1,N,ncmp);
                 for j=i+1:N
                     tmp = getSubspaceDistances(Spaces{i},Spaces{j},'projection');
                     Drow(1,j,:) = tmp;
                 end
                 DU(i,:,:) = Drow;
             end
             toc;
             D = DU;
             for i=1:1:ncmp
                D(:,:,i) = DU(:,:,i)+DU(:,:,i)';
             end
          % NP Manova test
          FStat = zeros(1,ncmp);
          for i=1:1:ncmp
               test = oneWayNPManova;
               test.D = D(:,:,i);
               test.n = [tb tb];
               test.t = 0;
               perform(test);
               FStat(i) = test.F;
           end
           if t<=0, return; end
           StatCount = false(ncmp,t);
           FStatpermC = zeros(t,length(FStat));
           nT = nX1+nX2;
           X = [X1; X2];
           disp('Permuting');
           tic;
           parfor i=1:t
                  ind = randperm(nT);
                  X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
                  X2perm = X(ind(nX1+1:end),:);
                  % Bootstrapping Subspaces
                  SpacesX1perm = cell(1,tb);
                  SpacesX2perm = cell(1,tb);
                  for k=1:tb
                      ind = randsample((1:nX1),nX1,true);
                      tmp = createSubspaceSVD(X1perm(ind,:)); %#ok<*PFBNS>
                      SpacesX1perm{k} = tmp(:,1:ncmp);
                      ind = randsample((1:nX2),nX2,true);
                      tmp = createSubspaceSVD(X2perm(ind,:));
                      SpacesX2perm{k} = tmp(:,1:ncmp);
                  end
                  % Generating the DistanceMatrix
                  Spacesperm = [SpacesX1perm';SpacesX2perm'];
                  N = size(Spacesperm,1);
                  DU = zeros(N,N,ncmp);
                  for k=1:N-1
                    Drow = zeros(1,N,ncmp);
                    for j=k+1:N
                        tmp = getSubspaceDistances(Spacesperm{k},Spacesperm{j},'projection');
                        Drow(1,j,:) = tmp;
                    end
                    DU(k,:,:) = Drow;
                  end
                  D = DU;
                  for k=1:1:ncmp
                      D(:,:,k) = DU(:,:,k)+DU(:,:,k)';
                  end
                  tmp = zeros(1,ncmp);
                  FStatperm = zeros(1,ncmp);
                  for k=1:1:ncmp
                      testperm = oneWayNPManova;
                      testperm.D = D(:,:,k);
                      testperm.n = [tb tb];
                      testperm.t = 0;
                      getSST(testperm);
                      getSSW(testperm);
                      getF(testperm);
                      FStatperm(k) = testperm.F;
                      tmp(k) = testperm.F>=FStat(k);
                  end
                  FStatpermC(i,:) = FStatperm;
                  StatCount(:,i) = tmp(:);
           end
           toc;
           figure;
           hold on;
           for i=1:1:t
               plot(FStatpermC(i,:),'r-.');
               drawnow;   
           end
           plot(FStat,'b','linewidth',3);
           out.pperm = (sum(StatCount,2)/t)';
           out.FStat = FStat;
           out.FstatpermC = FStatpermC;
end
% function out = MixedMomentFStatisticBoot(X1,X2,tb,t)
%          if nargin < 3, t = 0; end
%          nX1 = size(X1,1);
%          nX2 = size(X2,1);
%          % Bootstrapping Subspaces
%            SpacesX1 = cell(1,tb);
%            SpacesX2 = cell(1,tb);
%            disp('Booting Subspaces');
%            tic;
%            parfor i=1:tb
%                ind = randsample((1:nX1),round(nX1/2),true);
%                SpacesX1{i} = createSubspace(X1(ind,:)'); %#ok<*PFBNS>
%                ind = randsample((1:nX2),round(nX2/2),true);
%                SpacesX2{i} = createSubspace(X2(ind,:)');       
%            end
%            toc;
%          % Generating the DistanceMatrix
%              Spaces = [SpacesX1';SpacesX2'];
%              N = size(Spaces,1);
%              DU = zeros(N,N);
%              DL = zeros(N,N);
%              disp('Computing Distance Matrix');
%              tic;
%              parfor i=1:N-1
%                  Drow = zeros(1,N);
%                  for j=i+1:N
%                     Drow(j) = getSubspaceDistance(Spaces{i},Spaces{j},'projection frobenius');
%                  end
%                  DU(i,:) = Drow;
%                  DL(:,i) = Drow';
%              end
%              toc;
%              D = DU+DL;
%           % NP Manova test
%            test = oneWayNPManova;
%            test.D = D;
%            test.n = [tb tb];
%            test.t = t;
%            perform(test);
%            out.SpacesX1 = SpacesX1;
%            out.SpacesX2 = SpacesX2;
%            out.test = test;
%            out.FStat = test.F;
%            FStat = test.F;
%            if t<=0, return; end
%            StatCount = false(length(FStat),t);
%            nT = nX1+nX2;
%            X = [X1; X2];
%            disp('Permuting');
%            tic;
%            parfor i=1:t
%                   ind = randperm(nT);
%                   X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
%                   X2perm = X(ind(nX1+1:end),:);
%                   % Bootstrapping Subspaces
%                   SpacesX1perm = cell(1,tb);
%                   SpacesX2perm = cell(1,tb);
%                   for k=1:tb
%                       ind = randsample((1:nX1),round(nX1/2),true);
%                       SpacesX1perm{k} = createSubspace(X1perm(ind,:)'); %#ok<*PFBNS>
%                       ind = randsample((1:nX2),round(nX2/2),true);
%                       SpacesX2perm{k} = createSubspace(X2perm(ind,:)');       
%                   end
%                   % Generating the DistanceMatrix
%                   Spacesperm = [SpacesX1perm';SpacesX2perm'];
%                   N = size(Spacesperm,1);
%                   DU = zeros(N,N);
%                   DL = zeros(N,N);
%                   for k=1:N-1
%                     Drow = zeros(1,N);
%                     for j=k+1:N
%                         Drow(j) = getSubspaceDistance(Spacesperm{k},Spacesperm{j},'projection frobenius');
%                     end
%                     DU(k,:) = Drow;
%                     DL(:,k) = Drow';
%                   end
%                   D = DU+DL; 
%                   testperm = oneWayNPManova;
%                   testperm.D = D;
%                   testperm.n = [tb tb];
%                   testperm.t = 0;
%                   getSST(testperm);
%                   getSSW(testperm);
%                   getF(testperm);
%                   StatCount(i) = testperm.F>=FStat;
%            end
%            toc;
%            out.pperm = (sum(StatCount,2)/t)';
% end

% function out = MixedMomentFStatisticBootv2(X1,X2,tb,t,ncmp)
%          if nargin < 4, t = 0; end
%          if nargin < 5, ncmp = +inf; end
%          nX1 = size(X1,1);
%          nX2 = size(X2,1);
%          ncmp = min([nX1 nX2 ncmp]);
%          type = 'procrustes';
%          % Bootstrapping Subspaces
%            SpacesX1 = cell(1,tb);
%            SpacesX2 = cell(1,tb);
%            disp('Booting Subspaces');
%            tic;
%            parfor i=1:tb
%                ind = randsample((1:nX1),round(nX1/2),true);
%                tmp = createSubspace(X1(ind,:)'); %#ok<*PFBNS>
%                SpacesX1{i} = tmp(:,1:ncmp);
%                ind = randsample((1:nX2),round(nX2/2),true);
%                tmp = createSubspace(X2(ind,:)');       
%                SpacesX2{i} = tmp(:,1:ncmp);
%            end
%            toc;
%          % Generating the DistanceMatrix
%              Spaces = [SpacesX1';SpacesX2'];
%              N = size(Spaces,1);
%              DU = zeros(N,N,ncmp);
%              disp('Computing Distance Matrix');
%              tic;
%              parfor i=1:N-1
%                  Drow = zeros(1,N,ncmp);
%                  for j=i+1:N
%                      tmp = getSubspaceDistances(Spaces{i},Spaces{j},type);
%                      Drow(1,j,:) = tmp;
%                  end
%                  DU(i,:,:) = Drow;
%              end
%              toc;
%              D = DU;
%              for i=1:1:ncmp
%                 D(:,:,i) = DU(:,:,i)+DU(:,:,i)';
%              end
%           % NP Manova test
%           FStat = zeros(1,ncmp);
%           for i=1:1:ncmp
%                test = oneWayNPManova;
%                test.D = D(:,:,i);
%                test.n = [tb tb];
%                test.t = 0;
%                perform(test);
%                FStat(i) = test.F;
%            end
%            if t<=0, return; end
%            StatCount = false(ncmp,t);
%            DimCount = zeros(1,t);
%            nT = nX1+nX2;
%            X = [X1; X2];
%            disp('Permuting');
%            tic;
%            parfor i=1:t
%                   ind = randperm(nT);
%                   X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
%                   X2perm = X(ind(nX1+1:end),:);
%                   % Bootstrapping Subspaces
%                   SpacesX1perm = cell(1,tb);
%                   SpacesX2perm = cell(1,tb);
%                   Dimperm = +inf;
%                   for k=1:tb
%                       ind = randsample((1:nX1),round(nX1/2),true);
%                       SpacesX1perm{k} = createSubspace(X1perm(ind,:)'); %#ok<*PFBNS>
%                       Dimperm = min(Dimperm,size(SpacesX1perm{k},2));
%                       ind = randsample((1:nX2),round(nX2/2),true);
%                       SpacesX2perm{k} = createSubspace(X2perm(ind,:)');
%                       Dimperm = min(Dimperm,size(SpacesX2perm{k},2));
%                   end
%                   DimCount(i) = Dimperm;
%                   % Generating the DistanceMatrix
%                   Spacesperm = [SpacesX1perm';SpacesX2perm'];
%                   N = size(Spacesperm,1);
%                   DU = zeros(N,N,Dimperm);
%                   for k=1:N-1
%                     Drow = zeros(1,N,Dimperm);
%                     for j=k+1:N
%                         tmp = getSubspaceDistances(Spacesperm{k},Spacesperm{j},type);
%                         Drow(1,j,:) = tmp(1:Dimperm);
%                     end
%                     DU(k,:,:) = Drow;
%                   end
%                   D = DU;
%                   for k=1:1:Dimperm
%                       D(:,:,k) = DU(:,:,k)+DU(:,:,k)';
%                   end
%                   tmp = zeros(1,dim);
%                   for k=1:1:min(dim,Dimperm)
%                       testperm = oneWayNPManova;
%                       testperm.D = D(:,:,k);
%                       testperm.n = [tb tb];
%                       testperm.t = 0;
%                       getSST(testperm);
%                       getSSW(testperm);
%                       getF(testperm);
%                       tmp(k) = testperm.F>=FStat(k);
%                   end
%                   StatCount(:,i) = tmp(:);
%            end
%            FDim = min([dim DimCount]);
%            out.FDim = FDim;
%            out.FStat = FStat(1:FDim);
%            StatCount = StatCount(1:FDim,:);
%            toc;
%            out.pperm = (sum(StatCount,2)/t)';
% end
% % function out = MixedMomentFStatisticBoot(X1,X2,tb,t)
% %          if nargin < 3, t = 0; end
% %          nX1 = size(X1,1);
% %          nX2 = size(X2,1);
% %          % Bootstrapping Subspaces
% %            SpacesX1 = cell(1,tb);
% %            SpacesX2 = cell(1,tb);
% %            disp('Booting Subspaces');
% %            tic;
% %            parfor i=1:tb
% %                ind = randsample((1:nX1),round(nX1/2),true);
% %                SpacesX1{i} = createSubspace(X1(ind,:)'); %#ok<*PFBNS>
% %                ind = randsample((1:nX2),round(nX2/2),true);
% %                SpacesX2{i} = createSubspace(X2(ind,:)');       
% %            end
% %            toc;
% %          % Generating the DistanceMatrix
% %              Spaces = [SpacesX1';SpacesX2'];
% %              N = size(Spaces,1);
% %              DU = zeros(N,N);
% %              DL = zeros(N,N);
% %              disp('Computing Distance Matrix');
% %              tic;
% %              parfor i=1:N-1
% %                  Drow = zeros(1,N);
% %                  for j=i+1:N
% %                     Drow(j) = getSubspaceDistance(Spaces{i},Spaces{j},'projection frobenius');
% %                  end
% %                  DU(i,:) = Drow;
% %                  DL(:,i) = Drow';
% %              end
% %              toc;
% %              D = DU+DL;
% %           % NP Manova test
% %            test = oneWayNPManova;
% %            test.D = D;
% %            test.n = [tb tb];
% %            test.t = t;
% %            perform(test);
% %            out.SpacesX1 = SpacesX1;
% %            out.SpacesX2 = SpacesX2;
% %            out.test = test;
% %            out.FStat = test.F;
% %            FStat = test.F;
% %            if t<=0, return; end
% %            StatCount = false(length(FStat),t);
% %            nT = nX1+nX2;
% %            X = [X1; X2];
% %            disp('Permuting');
% %            tic;
% %            parfor i=1:t
% %                   ind = randperm(nT);
% %                   X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
% %                   X2perm = X(ind(nX1+1:end),:);
% %                   % Bootstrapping Subspaces
% %                   SpacesX1perm = cell(1,tb);
% %                   SpacesX2perm = cell(1,tb);
% %                   for k=1:tb
% %                       ind = randsample((1:nX1),round(nX1/2),true);
% %                       SpacesX1perm{k} = createSubspace(X1perm(ind,:)'); %#ok<*PFBNS>
% %                       ind = randsample((1:nX2),round(nX2/2),true);
% %                       SpacesX2perm{k} = createSubspace(X2perm(ind,:)');       
% %                   end
% %                   % Generating the DistanceMatrix
% %                   Spacesperm = [SpacesX1perm';SpacesX2perm'];
% %                   N = size(Spacesperm,1);
% %                   DU = zeros(N,N);
% %                   DL = zeros(N,N);
% %                   for k=1:N-1
% %                     Drow = zeros(1,N);
% %                     for j=k+1:N
% %                         Drow(j) = getSubspaceDistance(Spacesperm{k},Spacesperm{j},'projection frobenius');
% %                     end
% %                     DU(k,:) = Drow;
% %                     DL(:,k) = Drow';
% %                   end
% %                   D = DU+DL; 
% %                   testperm = oneWayNPManova;
% %                   testperm.D = D;
% %                   testperm.n = [tb tb];
% %                   testperm.t = 0;
% %                   getSST(testperm);
% %                   getSSW(testperm);
% %                   getF(testperm);
% %                   StatCount(i) = testperm.F>=FStat;
% %            end
% %            toc;
% %            out.pperm = (sum(StatCount,2)/t)';
% % end