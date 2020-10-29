function [Est] = myregressLOO(A,B,Model,type,index)
         n = size(A,1);
         Est = zeros(1,n);
         Ind = (1:n);
         %type = 'Euclidean';
         parfor i=1:n
            TrInd = setdiff(Ind,i);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            TestB = B(i,:);
            [TrA,TrB] = eliminateNAN(TrA,TrB);
            [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
            [~,~,Est(i),~] = getDistance(Model,[],TestB,type,TrM(1+index,:));
         end
end