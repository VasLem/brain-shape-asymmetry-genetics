function out = ShapeModelAnova2(X1,X2,AM,NL,NLsets,percvar,type)
         if nargin < 7, type = 'Mahalanobis'; end
         if nargin < 6, percvar = 99; end
         if nargin < 5, NLsets = 2; end
         if nargin < 4, NL = 0; end
         [n,nrV] = size(X1);        
         %Injection noise
         disp('Injecting noise');
         Set1 = zeros(n*NLsets,nrV);
         Set2 = zeros(n*NLsets,nrV);
         for i=1:n
             Set1((i-1)*3+1,:) = X1(i,:);
             Set2((i-1)*3+1,:) = X2(i,:);
             for k=1:1:NLsets-1
                 Set1((i-1)*3+1+k,:) = X1(i,:)+NL*randn(1,length(X1(i,:)));
                 Set2((i-1)*3+1+k,:) = X2(i,:)+NL*randn(1,length(X2(i,:)));
             end
         end
         % Building ShapeModel
         disp('Building Shape Model');
         SM = shapePCA;
         SM.RefScan = AM;
         getAverage(SM,[Set1',Set2']);
         getModel(SM,[Set1',Set2']);
         stripPercVar(SM,percvar);
         CSet1 = SM.Tcoeff(1:size(Set1,1),:);
         CSet2 = SM.Tcoeff(size(Set1,1)+1:end,:);
         if strcmp(type,'Mahalanobis')
             CSet1 = CSet1./repmat(SM.EigVal',n*NLsets,1);
             CSet2 = CSet2./repmat(SM.EigVal',n*NLsets,1);
         end
         out.ShapeModel = SM;
         out.CX1 = CSet1;
         out.CX2 = CSet2;
         [~,nrV] = size(CSet1);
         SS = zeros(4,nrV);
         tic;
         disp('Computing Anova');
         parfor i=1:nrV
            X = [CSet1(:,i) CSet2(:,i)];
            [~,TABLE] = anova2(X,NLsets,'off');
            SS(:,i) =  cell2mat(TABLE(2:5,2));
         end
         disp('Done');
         toc;
         % Individuals
         PC_I = SS(2,:);
         Total_I = sum(PC_I);
         PC_I = PC_I./((n-1));
         Total_I = Total_I/((n-1)*(nrV));        
         out.PC_I = PC_I;
         out.Total_I = Total_I;
         % Directional
         PC_D = SS(1,:);
         Total_D = sum(PC_D);
         Total_D = Total_D/(n-1);    
         out.PC_D = PC_D;
         out.Total_D = Total_D;
         % Fluctuating
         PC_F = SS(3,:);
         Total_F = sum(PC_F);
         PC_F = PC_F./((n-1));
         Total_F = Total_F/((n-1)*(nrV));       
         out.PC_F = PC_F;
         out.Total_F = Total_F;       
         % Error
         PC_E = SS(4,:);
         Total_E = sum(PC_E);
         PC_E = PC_E./((NLsets-1)*n*2);
         Total_E = Total_E/(((NLsets-1)*n*2*(nrV)));         
         out.PC_E = PC_E;
         out.Total_E = Total_E;
         
end

