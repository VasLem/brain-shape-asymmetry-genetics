function [NIndVar,Iter] = nestedBRIM(IndVar,DepVar,Model,type,index,nrIter)     
         n = size(IndVar,1);
         Ind = (1:n);
         nrBoot = length(index);
         if ~(size(IndVar,1)==size(DepVar,1)), error('different amount of observations'); end
         %Iter = zeros(size(IndVar,1),size(IndVar,2),nrIter+1);
         %Iter(:,:,1) = IndVar;
         %f = statusbar('nestedBRIM');drawnow;
         parfor i=1:n
            Iterfor = zeros(1,size(IndVar,2),nrIter+1);
            TrInd = setdiff(Ind,i);
            TrIndVar = IndVar(TrInd,:);
            TrDepVar = DepVar(TrInd,:);
            [TrIndVar,TrDepVar] = eliminateNAN(TrIndVar,TrDepVar);
            TestIndVar = IndVar(i,:);
            TestDepVar = DepVar(i,:);
            % Bootstrapping the training set
            TrIndVarNew = TrIndVar;
            for bootstep = 1:1:nrIter
                Iterfor(:,:,bootstep) = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,Model,type,index);
                for k=1:1:nrBoot
                    TrIndVarNew(:,index(k)) = myregressLOO(TrIndVar,TrDepVar,Model,type,index(k))';
                end 
                TrIndVar = TrIndVarNew;
            end
            Iterfor(:,:,end) = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,Model,type,index);
            Iter(i,:,:) = Iterfor;
            %statusbar(i/n,f);drawnow;
         end
         %delete(f);
         
         if ndims(Iter)==3
            NIndVar = Iter(:,:,end);
         else
            NIndVar = Iter(:,end);
         end
end

function [TestIndVarNew] = updateTestRIP(TrIndVar,TrDepVar,TestIndVar,TestDepVar,Model,type,index)
          TestIndVarNew = TestIndVar;
          [~,~,~,~,TrM,~,~,~] = plsregress(TrIndVar,TrDepVar,size(TrIndVar,2));
          for i=1:1:length(index)
              [~,~,TestIndVarNew(:,index(i)),~] = getDistance(Model,[],TestDepVar,type,TrM(1+index(i),:));          
          end
end