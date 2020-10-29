function [obj,NIndVar,Iter] = BRIM(obj,index,Cindex)
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
         % building correlationcoefficients
         [Xcor,Ycor] = eliminateNAN(obj.X,obj.Y);
         CorTemp = corrcoef([Xcor Ycor]);
         Cor = CorTemp(1:size(IndVar,2), size(IndVar,2)+1:size(CorTemp,2));
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
                Iterfor(:,:,bootstep) = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,obj,BootInd,Cor);
                TrIndVarNew(:,BootInd) = innerFoldBRIM(TrIndVar,TrDepVar,obj,BootInd,Cor);
                TrIndVar = TrIndVarNew;
             end
             Iterfor(:,:,end) = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,obj,BootInd,Cor);
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

function [TestIndVarNew] = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,obj,index,Cor)
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
%                         [~,~,TestIndVarNew(t,index(i))] = ...
%                         getDistance(obj.Model,[],TestDepVar(t,:),'mahalanobis',TrM(1+index(i),:));
                          %disp('Ok');
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
                    case 'PLSR Weighted'
                     for t=1:1:nrTest  
                       for i=1:1:length(index) 
                        [~,~,TestIndVarNew(t,index(i))] = ...
                        getDistanceWeighted(obj.Model,[],TestDepVar(t,:),'weighted',TrM(1+index(i),:),Cor(i,:));
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

function [Est] = innerFoldBRIM(A,B,obj,index,Cor)
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
                    case 'PLSR Weighted'
                        for j=1:1:nrTest
                            TestB = B(TestInd(j),:);
                            for k=1:1:nr2Boot
                                [~,~,Est(TestInd(j),k),~] = getDistanceWeighted(obj.Model,[],TestB,'weighted',TrM(1+index(k),:),Cor(k,:));
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

% function out = getDistancev2(obj,A,type,Direction)
%          B = repmat(obj.Model.AvgCoeff',size(A,1),1);% to be changed by AvgY
%          C = repmat(Direction,size(A,1),1);
%          D = A-B;
%          switch lower(type)
%              case 'euclidean'
%                 out = sqrt(sum((D).^2,2));
%              case 'mahalanobis'
%                 out = sqrt(sum(((D)./repmat(obj.Model.EigStd',size(A,1),1)).^2,2));
%          end
%          normD = sqrt(sum((D.^2),2));
%          D = D./repmat(normD,1,size(D,2));
%          switch lower(type)
%              case 'euclidean'
%                 % do nothing
%              case 'mahalanobis'
%                 D = D./repmat(obj.Model.EigStd',size(D,1),1);
%                 C = C./repmat(obj.Model.EigStd',size(C,1),1);
%          end
%          T = diag(C*D');
%          P1 = diag(C*C');P2 = diag(D*D');
%          N = sqrt(P1.*P2);
%          angles = T./N;
%          out = angles.*out;
% end