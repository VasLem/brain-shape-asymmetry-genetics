function out = FirstOrderMomentFStatisticPaired(X1,X2,t)
         if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting Centroids
           AvgX1 = mean(X1);
           AvgX2 = mean(X2);
         % generating Effect output  
           out.AvgX1 = AvgX1;
           out.AvgX2 = AvgX2;
           out.Difference = AvgX1-AvgX2;
           out.Distance = sqrt(sum(out.Difference.^2));
         % generating DistanceMatrix + NP Manova test
           test = oneWayNPManova;
           test.D = squareform(pdist([X1;X2],'euclidean'));
           test.n = [nX1 nX2];
           test.t = 0;
           perform(test);
           FStat = test.F;
           out.FStat = FStat;
           if t<=0, return; end
           % generating test-statistic (Euclidean Distance)
           StatCount = false(1,t);
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
                  testperm = oneWayNPManova;
                  testperm.D = squareform(pdist([X1perm;X2perm],'euclidean'));
                  testperm.n = [nX1 nX2];
                  testperm.t = 0;
                  getSST(testperm);
                  getSSW(testperm);
                  getF(testperm);
                  StatCount(i) = testperm.F>=FStat;
           end
           toc;
           out.pperm = sum(StatCount)/t;     
end