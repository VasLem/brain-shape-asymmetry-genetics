function [TestIndVarNew] = updateMultipleTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,Model,type,index)
          [~,~,~,~,TrM,~,~,~] = plsregress(TrIndVar,TrDepVar,size(TrIndVar,2));
          nrTest = size(TestIndVar,1);
          TestIndVarNew = TestIndVar;
          for t=1:1:nrTest
            for i=1:1:length(index)
                [~,~,TestIndVarNew(t,index(i))] = getDistance(Model,[],TestDepVar(t,:),type,TrM(1+index(i),:));          
            end
          end
end