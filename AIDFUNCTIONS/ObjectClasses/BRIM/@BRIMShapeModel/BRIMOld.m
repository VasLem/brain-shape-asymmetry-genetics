function [obj,NIndVar,Iter] = BRIMOld(obj,index,Cindex)
         if nargin < 3, Cindex = [];end
         if ~checkNrObservations(obj), error('different amount of observations between X and Y'); end
         % Building Independent and Dependent Variables
         n = obj.nrO;
         Ind = (1:n);
         nr2Condition = length(Cindex);
         nr2Boot = length(index);
         if ~nr2Condition==0
            IndVar = [obj.RIPX(:,Cindex) obj.X(:,index)];
         else
            IndVar = obj.X(:,index);
         end       
         BootInd = (nr2Condition+1:nr2Condition+nr2Boot);
         DepVar = obj.Y;
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
%                 TrIndVarNew(:,BootInd) = innerFoldBRIM(TrIndVar,TrDepVar,obj,BootInd);
                for k=1:1:nr2Boot% inner fold
                    TrIndVarNew(:,BootInd(k)) = innerFoldBRIM(TrIndVar,TrDepVar,obj,BootInd(k))';
                end% end inner fold 
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
         % linking results to BRIMShapeModel
         obj.RIPX(:,index) = NIndVar(:,BootInd);
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
                        getDistance(obj.Model,[],TestDepVar(t,:),'Mahalanobis',TrM(1+index(i),:));          
                       end
                      end
                    case 'PLSR Euclidean'
                     for t=1:1:nrTest  
                       for i=1:1:length(index) 
                        [~,~,TestIndVarNew(t,index(i))] = ...
                        getDistance(obj.Model,[],TestDepVar(t,:),'Euclidean',TrM(1+index(i),:));
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
         Est = zeros(1,n);
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
                            [~,~,Est(TestInd(j)),~] = getDistance(obj.Model,[],TestB,'Mahalanobis',TrM(1+index,:));
                        end
                    case 'PLSR Euclidean'
                        for j=1:1:nrTest
                            TestB = B(TestInd(j),:);
                            [~,~,Est(TestInd(j)),~] = getDistance(obj.Model,[],TestB,'Euclidean',TrM(1+index,:));
                        end
                    case 'PLSR PINV'
%                         TrMinv = pinv(TrM(2:end,:));
%                         avg = mean(TrA);
%                         tmp = B(TestInd,:)*TrMinv+repmat(avg,nrTest,1);
                    otherwise
            end
         end% end Inner Fold
end