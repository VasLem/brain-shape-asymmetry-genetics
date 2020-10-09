function out = ShapeModelAnovaN(X,G,AM,percvar,type,varargin)
         [n,nrV] = size(X);
         % Building ShapeModel
         tic;
         disp('Building Shape Model');
         SM = shapePCA;
         SM.RefScan = AM;
         getAverage(SM,X');
         getModel(SM,X');
         stripPercVar(SM,percvar);
         if strcmpi(type,'mahalanobis')
             X = SM.Tcoeff./repmat(SM.EigVal',n,1);
         else
             X = SM.Tcoeff;
         end
         [~,nrV] = size(X);
         % dummy test to determine sizes for allocation
         [~, TABLE] = anovan(X(:,1), G, varargin{:},'display','off');
         nE = size(TABLE,1)-2;
         DF = cell2mat(TABLE(2:end-1,3));
         Names = TABLE(2:end-1,1);
         SS = zeros(nE,nrV);
         SS(:,1) =  cell2mat(TABLE(2:end-1,2));
         parfor i=2:nrV
            [~,TABLE] = anovan(X(:,i),G,varargin{:},'display','off');
            SS(:,i) =  cell2mat(TABLE(2:end-1,2));
         end
         for i=1:1:nE
             Effect(i).Name = Names{i}; %#ok<*AGROW>
             Effect(i).PC = SS(i,:);
%              Effect(i).PC = sum(Effect(i).SS);
             Effect(i).Total = sum(Effect(i).PC);
             Effect(i).PCD = sqrt(Effect(i).PC)/DF(i);
             Effect(i).PC = Effect(i).PC./(DF(i));
             Effect(i).PCD = sqrt(Effect(i).PC);
             Effect(i).Total = Effect(i).Total/((nrV)*DF(i));           
         end
         out.Effect = Effect;
         out.SM = SM;
         toc;
end