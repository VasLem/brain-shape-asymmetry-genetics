function [F,P] = FStatTest(X1,X2,X3,t)
         % initializing  
           nX1 = size(X1,1);nX2 = size(X2,1);nX3 = size(X3,1);nX = nX1+nX2+nX3;
           X = [X1;X2;X3];
         % group averages
           AvgX1 = mean(X1);AvgX2 = mean(X2);AvgX3 = mean(X3);AvgX = mean(X);   
         % SSB
           SSB = sum(sum(([AvgX1;AvgX2;AvgX3] - repmat(AvgX,3,1)).^2,2).*[nX1;nX2;nX3])/2;
         % SSW
           SSW = sum([sum((X1-repmat(AvgX1,nX1,1)).^2,2);sum((X2-repmat(AvgX2,nX2,1)).^2,2);sum((X3-repmat(AvgX3,nX3,1)).^2,2)])/(nX-2);
         % F
           F = SSB/SSW;
         % Permutation test  
           StatCount = false(1,t);
           for i=1:t
                 % initialize permutation
                  ind = randperm(nX);
                  X1perm = X(ind(1:nX1),:);X2perm = X(ind(nX1+1:nX1+nX2),:);X3perm = X(ind(nX1+nX2+1:end),:); %#ok<*PFBNS>
                 % group averages 
                  AvgX1perm = mean(X1perm);AvgX2perm = mean(X2perm);AvgX3perm = mean(X3perm);AvgXperm = mean(X(ind,:));
                 % SSB
                  SSBperm = sum(sum(([AvgX1perm;AvgX2perm;AvgX3perm] - repmat(AvgXperm,3,1)).^2,2).*[nX1;nX2;nX3])/2;
                 % SSW 
                  SSWperm = sum([sum((X1perm-repmat(AvgX1perm,nX1,1)).^2,2);sum((X2perm-repmat(AvgX2perm,nX2,1)).^2,2);sum((X3perm-repmat(AvgX3perm,nX3,1)).^2,2)])/(nX-2);
                 % F
                  StatCount(i) =  (SSBperm/SSWperm)>=F;
           end
           P = sum(StatCount)/t;
end