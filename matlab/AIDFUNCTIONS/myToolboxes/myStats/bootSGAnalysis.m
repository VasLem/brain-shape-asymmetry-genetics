function [out] = bootSGAnalysis(A,B,Model,t)
         if nargin < 4, t = 100; end
         n = size(A,1);
         sampleauc = zeros(1,t);
         samplex = zeros(n,t);
         sampley = zeros(n,t);
         Ind = (1:n);
         parfor b=1:t
            ind = randsample((1:n),n,true);
            Afor = A(ind,:);
            Bfor = B(ind,:);
            Estfor = zeros(1,n);
            for i=1:1:n       
                TrInd = setdiff(Ind,i);
                TrA = Afor(TrInd,:);
                TrB = Bfor(TrInd,:);
                TestB = Bfor(i,:);
                [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
                [~,~,Estfor(i),~] = getDistance(Model,[],TestB,'Mahalanobis',TrM(end,:));
            end
            X = Afor(:,end);Y = Estfor;
            X(find(X == -1)) = 0;
            [samplex(:,b),sampley(:,b),sampleauc(b)] = roc_calc(~logical(X),Y);
         end
         sampleauc = sort(sampleauc);
         out.upper = sampleauc(0.975*t);
         out.lower= sampleauc(0.025*t+1);
         out.avg = mean(sampleauc);
         out.sampleAUC = sampleauc;
         out.samplex = samplex;
         out.sampley = sampley;
end


% parfor b=1:t
%             ind = randsample((1:n),n,true);
%             Afor = A(ind,:);
%             Bfor = B(ind,:);
%             Estfor = zeros(1,n);
%             for i=1:1:n
%                 TrInd = setdiff(Ind,i);
%                 TrA = Afor(TrInd,:);
%                 TrB = Bfor(TrInd,:);
%                 TestB = Bfor(i,:);
%                 [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
%                 [~,~,Estfor(i),~] = getDistance(Model,[],TestB,'Mahalanobis',TrM(end,:));
%             end
%             X = Afor(:,end);Y = Estfor;
%             X(find(X == -1)) = 0;
%             [sx,sy,sample(b)] = roc_calc(~logical(X),Y);
%          end