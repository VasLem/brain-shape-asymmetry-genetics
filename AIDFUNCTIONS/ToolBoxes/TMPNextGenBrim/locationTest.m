function [D,P] = locationTest(X1,X2,t)
         % initializing  
           nX1 = size(X1,1);nX2 = size(X2,1);
         % Getting First Order moments & Residus
           AvgX1 = mean(X1);AvgX2 = mean(X2);
         % Getting Distance Based Test Statistic
           D =  sqrt((AvgX1-AvgX2)*(AvgX1-AvgX2)');
         % Permutation test  
           StatCount = false(1,t);
           nT = nX1+nX2;
           X = [X1; X2];
           for i=1:t
                  ind = randperm(nT);
                  X1perm = X(ind(1:nX1),:);X2perm = X(ind(nX1+1:end),:); %#ok<*PFBNS>
                  AvgX1perm = mean(X1perm);AvgX2perm = mean(X2perm);
                  Dperm = sqrt((AvgX1perm-AvgX2perm)*(AvgX1perm-AvgX2perm)');
                  StatCount(i) = Dperm>=D;
           end
           P = sum(StatCount)/t;           
end