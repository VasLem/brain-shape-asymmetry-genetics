function out = FirstOrderMomentDStatisticPaired(X1,X2,t)
         if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting First Order moments & Residus
           AvgX1 = mean(X1);
           AvgX2 = mean(X2);
           ResX1 = X1-repmat(AvgX1,nX1,1);
           ResX2 = X2-repmat(AvgX2,nX2,1);
         % Getting Distance Based Test Statistic
           DStat =  sqrt((AvgX1-AvgX2)*(AvgX1-AvgX2)');
         % Getting within Sample Distance standard deviations 
           s1 = sqrt(sum((sum(ResX1.^2,2)))/(nX1-1));
           s2 = sqrt(sum((sum(ResX2.^2,2)))/(nX2-1));
         % Getting pooled standard deviation
           s = sqrt(((nX1-1)*s1^2+(nX2-1)*s2^2)/(nX1+nX2));
         % Getting Cohen's Distance
           CohenD = DStat/s;
         % generating Effect output  
           out.AvgX1 = AvgX1;
           out.AvgX2 = AvgX2;
           out.Difference = AvgX1-AvgX2;
           out.Dstat = DStat;
           out.s1 = s1;
           out.s2 = s2;
           out.s = s;
           out.CohenD = CohenD;
         % Permutation test  
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
                  AvgX1perm = mean(X1perm);
                  AvgX2perm = mean(X2perm);
                  Dperm = sqrt((AvgX1perm-AvgX2perm)*(AvgX1perm-AvgX2perm)');
                  StatCount(i) = Dperm>=DStat;
           end
           toc;
           out.pperm = sum(StatCount)/t;           
end