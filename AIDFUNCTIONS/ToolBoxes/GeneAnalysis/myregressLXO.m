function [Est] = myregressLXO(A,B,Model,type,index,X)
         n = size(A,1);
         Est = zeros(1,n);
         Ind = (1:n);
         if X==0, X=n; end
         K = round(n/X);
         Indices = crossvalind('Kfold',n,K);
         for i=1:K
            TestInd = find(Indices==i);
            nrTest = length(TestInd);
            TrInd = setdiff(Ind,TestInd);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            [TrA,TrB] = eliminateNAN(TrA,TrB);
            [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
            for j=1:1:nrTest
                TestB = B(TestInd(j),:);
                [~,~,Est(TestInd(j)),~] = getDistance(Model,[],TestB,type,TrM(1+index,:));
            end
         end
end