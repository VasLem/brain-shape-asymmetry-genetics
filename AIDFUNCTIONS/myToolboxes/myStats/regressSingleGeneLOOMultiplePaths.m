function [Est,CorEst,pCorEst] = regressSingleGeneLOOMultiplePaths(A,B,Model,type)      
         % additive path
         %disp('Additive');
         index = (1:size(A,1));
         Est(1,:) = myregressSingleGeneLOO(A,B,Model,type,index);
         % Homozygotes only path
         %disp('Homozygotes only');
         index = find(A(:,3));
         Est(2,:) = myregressSingleGeneLOO(A,B,Model,type,index);
         % added to minus one
         %disp('Added minus one');
         index = (1:size(A,1));
         Atemp = A(:,3);
         Atemp(find(A(:,3)==0)) = -1;
         Atemp = [A(:,2) Atemp];
         Est(3,:) = myregressSingleGeneLOO(Atemp,B,Model,type,index);
         % added to plus one
         %disp('Added plus one');
         index = (1:size(A,1));
         Atemp = A(:,3);
         Atemp(find(A(:,3)==0)) = 1;
         Atemp = [A(:,2) Atemp];
         Est(4,:) = myregressSingleGeneLOO(Atemp,B,Model,type,index);
         Est(5,:) = A(:,3)';
         % Correlations
         if nargout == 1, return; end
         CorEst = zeros(5,5);
         pCorEst = zeros(5,5);
         for i=1:1:5
            for j = 1:1:5
                [c,p] = corrcoef(Est(i,:)',Est(j,:)');
                CorEst(i,j) = c(1,2);
                pCorEst(i,j) = p(1,2);
            end            
         end
end

function [Est] = myregressSingleGeneLOO(A,B,Model,type,index)
         n = size(A,1);
         Est = zeros(1,n);
         Ind = (1:n);
         %type = 'Euclidean';
         for i=1:1:n
            TrInd = setdiff(Ind,i);
            TrInd = intersect(TrInd,index);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            TestB = B(i,:);
            [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
            [~,~,Est(i),~] = getDistance(Model,[],TestB,type,TrM(end,:));
         end
end