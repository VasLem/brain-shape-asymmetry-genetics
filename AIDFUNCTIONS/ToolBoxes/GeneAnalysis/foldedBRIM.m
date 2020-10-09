function [NIndVar,Iter] = foldedBRIM(IndVar,DepVar,Model,type,index,nrIter,X,Xi)
         if nargin < 8, Xi = X; end
         n = size(IndVar,1);
         Ind = (1:n);
         nr2Boot = length(index);
         if ~(size(IndVar,1)==size(DepVar,1)), error('different amount of observations'); end
         % subdivide the dataset into subsets
         K = round(n/X);
         Indices = crossvalind('Kfold',n,K);
         ForResults = cell(1,K);
         parfor i=1:K
            TestInd = find(Indices==i);
            nrTest = length(TestInd);
            TrInd = setdiff(Ind,TestInd);          
            %TrInd = setdiff(Ind,i);
            TrIndVar = IndVar(TrInd,:);
            TrDepVar = DepVar(TrInd,:);
            [TrIndVar,TrDepVar] = eliminateNAN(TrIndVar,TrDepVar);
            %TestIndVar = IndVar(i,:);
            %TestDepVar = DepVar(i,:);
            TestIndVar = IndVar(TestInd,:);
            TestDepVar = DepVar(TestInd,:);
            Iterfor = zeros(nrTest,size(IndVar,2),nrIter+1);
            % Bootstrapping the training set
            TrIndVarNew = TrIndVar;
            for bootstep = 1:1:nrIter
                Iterfor(:,:,bootstep) = updateMultipleTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,Model,type,index);
                for k=1:1:nr2Boot
                    TrIndVarNew(:,index(k)) = myregressLXO(TrIndVar,TrDepVar,Model,type,index(k),Xi)';
                end 
                TrIndVar = TrIndVarNew;
            end
            Iterfor(:,:,end) = updateMultipleTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,Model,type,index);
            ForResults{i} = Iterfor;
         end
         % assemble results
         Iter = zeros(size(IndVar,1),size(IndVar,2),nrIter+1);
         for i=1:1:K
             TestInd = find(Indices==i);
             Iter(TestInd,:,:) = ForResults{i};
         end
         Iter = squeeze(Iter);
         if ndims(Iter)==3
            NIndVar = Iter(:,:,end);
         else
            NIndVar = Iter(:,end);
         end
end

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