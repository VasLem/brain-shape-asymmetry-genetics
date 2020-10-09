function Effect = ProcrustesAnovaN(X,G,varargin)
         [~,nrV] = size(X);
         % dummy test to determine sizes for allocation
         [~, TABLE] = anovan(X(:,1), G, varargin{:},'display','off');
         nE = size(TABLE,1)-2;
         DF = cell2mat(TABLE(2:end-1,3));
         Names = TABLE(2:end-1,1);
         SS = zeros(nE,nrV);
         SS(:,1) =  cell2mat(TABLE(2:end-1,2));
         tic;
         parfor i=2:nrV
            [~,TABLE] = anovan(X(:,i),G,varargin{:},'display','off'); %#ok<PFBNS>
            SS(:,i) =  cell2mat(TABLE(2:end-1,2));
         end
         toc;
         for i=1:1:nE
             Effect(i).Name = Names{i}; %#ok<*AGROW>
             Effect(i).SS = reshape(SS(i,:),3,length(SS(i,:))/3);
             Effect(i).LM = sum(Effect(i).SS);
             Effect(i).Total = sum(Effect(i).LM);
             Effect(i).LMD = sqrt(Effect(i).LM/DF(i));
             Effect(i).LM = Effect(i).LM./(3*DF(i));
             %Effect(i).LMD = sqrt(Effect(i).LM);
             Effect(i).Total = Effect(i).Total/((nrV-7)*DF(i));
             Effect(i).df = DF(i);
         end      
end