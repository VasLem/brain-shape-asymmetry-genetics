function out = classifierInfoGain(EstClass,TrueClass)
         EstClass = EstClass(:)';
         TrueClass = TrueClass(:)';
         index = find(~isnan(TrueClass));
         EstClass = EstClass(index);
         TrueClass = TrueClass(index);
         CatLabel = unique(TrueClass);
         nrCat = length(CatLabel);
         nrO = length(TrueClass);
         CatP = zeros(1,nrO);
         E = 0;
         for i=1:1:nrCat
             P = sum(TrueClass==CatLabel(i))/length(TrueClass);
             E = E + P*log2(P); 
             CatP(TrueClass==CatLabel(i)) = P;
         end
         E = -1*E;
         res = (EstClass==TrueClass);
         correct = find(res==1);
         nrCorrect = length(correct);
         miss = find(res==0);
         nrMiss = length(miss);              
         Ic =  -1*log2(CatP(correct)) + log2(ones(1,nrCorrect));
         Im = -1*(-1*log2(1-CatP(miss))+log2(ones(1,nrMiss)));     
         Ia = sum([Ic Im])/nrO;
         Ir = Ia/E;
         out.Ia = Ia;
         out.Ir = Ir;
         out.E = E;
end