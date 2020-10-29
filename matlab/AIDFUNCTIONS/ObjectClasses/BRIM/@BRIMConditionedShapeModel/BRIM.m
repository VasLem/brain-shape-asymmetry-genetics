function [NIndVar,Iter] = BRIM(obj,IndVar,DepVar,index,Cindex)
         if nargin < 3, Cindex = [];end
         if ~(size(IndVar,1)==size(DepVar,1)), error('different amount of observations'); end
         % Building Independent and Dependent Variables
         n = size(IndVar,1);
         Ind = (1:n);
         nr2Condition = length(Cindex);
         nr2Boot = length(index);
         BootInd = (nr2Condition+1:nr2Condition+nr2Boot);
         % subdivide the dataset into subsets
         K = round(n/obj.TrainingOuterFold);
         Indices = crossvalind('Kfold',n,K);
         ForResults = cell(1,K);
         parfor i=1:K% Outer Fold
             TestInd = find(Indices==i);
             nrTest = length(TestInd);
             TrInd = setdiff(Ind,TestInd);
             TrIndVar = IndVar(TrInd,:); %#ok<*PFBNS>
             TrDepVar = DepVar(TrInd,:);
             [TrIndVar,TrDepVar] = eliminateNAN(TrIndVar,TrDepVar);
             TestIndVar = IndVar(TestInd,:);
             TestDepVar = DepVar(TestInd,:);
             Iterfor = zeros(nrTest,size(IndVar,2),obj.TrainingIterations+1);
             % Bootstrapping the training set
             TrIndVarNew = TrIndVar;
             for bootstep = 1:1:obj.TrainingIterations
                Iterfor(:,:,bootstep) = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,obj,BootInd);
                TrIndVarNew(:,BootInd) = innerFoldBRIM(TrIndVar,TrDepVar,obj,BootInd);
                TrIndVar = TrIndVarNew;
             end
             Iterfor(:,:,end) = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,obj,BootInd);
             ForResults{i} = Iterfor;
         end% end Outer Fold
         % assemble results
         Iter = zeros(size(IndVar,1),size(IndVar,2),obj.TrainingIterations+1);
         for i=1:1:K
             TestInd = find(Indices==i);
             Iter(TestInd,:,:) = ForResults{i}; %#ok<*FNDSB>
         end
         Iter = squeeze(Iter);
         if ndims(Iter)==3
            NIndVar = Iter(:,:,end);
         else
            NIndVar = Iter(:,end);
         end
end

function [TestIndVarNew] = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,obj,index)
          switch obj.BRIMModelType
              case 'PLSR'
                [~,~,~,~,TrM,~,~,~] = plsregress(TrIndVar,TrDepVar,size(TrIndVar,2));
              otherwise
          end
          nrTest = size(TestIndVar,1);
          TestIndVarNew = TestIndVar;
          switch obj.RIPEstimationMethod
                    case 'PLSR Mahalanobis'
                      for t=1:1:nrTest  
                       for i=1:1:length(index)
                        [~,~,TestIndVarNew(t,index(i))] = ...
                        getDistance(obj.Model,[],TestDepVar(t,:),'mahalanobis',TrM(1+index(i),:));          
                       end
                      end
                    case 'PLSR Euclidean'
                     for t=1:1:nrTest  
                       for i=1:1:length(index) 
                        [~,~,TestIndVarNew(t,index(i))] = ...
                        getDistance(obj.Model,[],TestDepVar(t,:),'euclidean',TrM(1+index(i),:));
                       end
                     end
                    case 'PLSR PINV'
                        TrMinv = pinv(TrM(2:end,:));
                        avg = mean(TrIndVar);
                        tmp = TestDepVar*TrMinv+repmat(avg,nrTest,1);
                        TestIndVarNew(:,index) = tmp(:,index);
                    otherwise
          end
end

function [Est] = innerFoldBRIM(A,B,obj,index)
         n = size(A,1);
         nr2Boot = length(index);
         Est = zeros(n,nr2Boot);
         Ind = (1:n);
         K = round(n/obj.TrainingInnerFold);
         Indices = crossvalind('Kfold',n,K);
         for i=1:K% Inner Fold
            TestInd = find(Indices==i);
            nrTest = length(TestInd);
            TrInd = setdiff(Ind,TestInd);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            [TrA,TrB] = eliminateNAN(TrA,TrB);
            switch obj.BRIMModelType
                case 'PLSR'
                     [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
                otherwise
            end
            switch obj.RIPEstimationMethod
                    case 'PLSR Mahalanobis'
                        for j=1:1:nrTest
                            TestB = B(TestInd(j),:);
                            for k=1:1:nr2Boot
                                [~,~,Est(TestInd(j),k),~] = getDistance(obj.Model,[],TestB,'mahalanobis',TrM(1+index(k),:));
                            end
                        end
                    case 'PLSR Euclidean'
                        for j=1:1:nrTest
                            TestB = B(TestInd(j),:);
                            for k=1:1:nr2Boot
                                [~,~,Est(TestInd(j),k),~] = getDistance(obj.Model,[],TestB,'euclidean',TrM(1+index(k),:));
                            end
                        end
                    case 'PLSR PINV'
                        TrMinv = pinv(TrM(2:end,:));
                        avg = mean(TrA);
                        tmp = B(TestInd,:)*TrMinv+repmat(avg,nrTest,1);
                        Est(TestInd,:) = tmp(:,index);
                    otherwise
            end
         end% end Inner Fold
end