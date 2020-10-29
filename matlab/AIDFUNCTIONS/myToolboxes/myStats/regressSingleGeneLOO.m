function [Est] = regressSingleGeneLOO(A,B,Model,type)
         n = size(A,1);
         Est = zeros(1,n);
         Ind = (1:n);
         %type = 'Euclidean';
         for i=1:1:n
            TrInd = setdiff(Ind,i);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            TestB = B(i,:);
            [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
            [~,~,Est(i),~] = getDistance(Model,[],TestB,type,TrM(end,:));
         end
end