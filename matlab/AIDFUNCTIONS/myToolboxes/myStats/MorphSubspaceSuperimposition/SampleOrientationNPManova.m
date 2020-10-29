function out = SampleOrientationNPManova(X1,X2,tb,tp)
         if nargin < 3, tp = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting Subspace centroids & Residus
             AvgX1 = mean(X1);
             AvgX2 = mean(X2);
             ResX1 = X1;%-repmat(AvgX1,nX1,1);
             ResX2 = X2;%-repmat(AvgX2,nX2,1);
         % Getting Subspace Scale
             SX1 = mean(sqrt(sum(ResX1.^2,2)));
             SX2 = mean(sqrt(sum(ResX2.^2,2)));    
             %ResX1 = ResX1/SX1;
             %ResX2 = ResX2/SX2;
         % Bootstrapping Subspaces
             SpacesX1 = cell(1,tb);
             %SpacesX1{1} = orth(ResX1');
             SpacesX2 = cell(1,tb);
             %SpacesX2{1} = orth(ResX2');
             disp('Booting Subspaces');
             tic;
             parfor i=1:tb
                 ind = randsample((1:nX1),round(nX1/2),true);
                 SpacesX1{i} = createSubspace(ResX1(ind,:)'); %#ok<*PFBNS>
                 ind = randsample((1:nX2),round(nX2/2),true);
                 SpacesX2{i} = createSubspace(ResX2(ind,:)');       
             end
             toc;
         % Generating the DistanceMatrix
             Spaces = [SpacesX1';SpacesX2'];
             N = size(Spaces,1);
             DU = zeros(N,N);
             DL = zeros(N,N);
             disp('Computing Distance Matrix');
             tic;
             parfor i=1:N-1
                 Drow = zeros(1,N);
                 for j=i+1:N
                    Drow(j) = getSubspaceDistance(Spaces{i},Spaces{j},'projection frobenius');
                 end
                 DU(i,:) = Drow;
                 DL(:,i) = Drow';
             end
             toc;
             D = DU+DL;
         % NP Manova test
           test = oneWayNPManova;
           test.D = D;
           test.n = [tb tb];
           test.t = tp;
           perform(test);
           out.SpacesX1 = SpacesX1;
           out.SpacesX2 = SpacesX2;
           out.Test = test;
end



