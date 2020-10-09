function [out,figures] = permBootLOOSGNullANOVAROCAnalysis(A,B,Model,t,string,type,G)

%          A = AABfor;
%          B = BABfor;
%          t = 0;
%          string = name;
%          type = 'Mahalanobis';
         % Retrieving Estimates in a LOO scenario      
         n = size(A,1);
         Est = zeros(1,n);
         Ind = (1:n);
         parfor i=1:n
            TrInd = setdiff(Ind,i);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            TestB = B(i,:);
            TrM = myBootstrapRegressSingleGene(TrA,TrB,Model,type);
            [~,~,Est(i),~] = getDistance(Model,[],TestB,type,TrM);
         end
         X = G;Y = Est;
         out.GeneDP = Est;
         %out = [];figures = []; return;
         % perform ANOVA
         %[FA,pFA] = myAnova(X,Y);
         [out.ANOVA,figures.figB] = dGPANOVAAnalysis(X,Y,t,string);
%          out.ANOVA.F = FA;
%          out.ANOVA.pF = pFA;
         % Perform Homozygote ROC
         msH = 1;msM = 1;msP = 1;change = true;
         out.ROC = myROCtests(X,Y,msH,msM,msP,change);
         out.ROC.P.pauc = normpdf((out.ROC.P.auc-0.5)/out.ROC.P.se,0,1);
         out.ROC.M.pauc = normpdf((out.ROC.M.auc-0.5)/out.ROC.M.se,0,1);
         out.ROC.H.pauc = normpdf((out.ROC.H.auc-0.5)/out.ROC.H.se,0,1);
%          if t==0, figures = []; return; end
%          disp('Permuting');
%          n = size(A,1);
%          FAnull = zeros(1,t);
%          Haucnull = zeros(1,t);
%          Hxnull = zeros(length(out.ROC.H.x),t);
%          Hynull = zeros(length(out.ROC.H.y),t);
%          Maucnull = zeros(1,t);
%          Mxnull = zeros(n,t);
%          Mynull = zeros(n,t);
%          Paucnull = zeros(1,t);
%          Pxnull = zeros(n,t);
%          Pynull = zeros(n,t);
%          change = false;
%          parfor b=1:t
%             ind = randperm(n);
%             Afor = A;
%             Afor(:,end) = Afor(ind,end);
%             Gfor = G(ind);
%             Bfor = B;       
%             %Estfor = myregressSingleGeneLOO(Afor,Bfor,Model,type);
%             Estfor = myBootstrapRegressSingleGeneLOO(Afor,Bfor,Model,type);
%             Xfor = Gfor;Yfor = Estfor;
%             [FAnull(b),~] = myAnova(Xfor,Yfor);
%             
%             ROCfor = myROCtests(Xfor,Yfor,msH,msM,msP,change);
%             
%             Haucnull(b) = ROCfor.H.auc;
%             Hxnull(:,b) = ROCfor.H.x;
%             Hynull(:,b) = ROCfor.H.y;
%             
%             Maucnull(b) = ROCfor.M.auc;
%             Mxnull(:,b) = ROCfor.M.x;
%             Mynull(:,b) = ROCfor.M.y;
%             
%             Paucnull(b) = ROCfor.P.auc;
%             Pxnull(:,b) = ROCfor.P.x;
%             Pynull(:,b) = ROCfor.P.y;
%             
%          end
%          out.ANOVA.ppermF = sum(FAnull>=FA)/t;
%          out.ROC.H.ppermauc = sum(Haucnull>=out.ROC.H.auc)/t;
%          out.ROC.M.ppermauc = sum(Maucnull>=out.ROC.M.auc)/t;
%          out.ROC.P.ppermauc = sum(Paucnull>=out.ROC.P.auc)/t;
%          if nargout == 1, return; end
%          figures.figH = plotROCresults(out.ROC.H,Hxnull,Hynull,[string '_H '],t);
%          figures.figM = plotROCresults(out.ROC.M,Mxnull,Mynull,[string '_M '],t);
%          figures.figP = plotROCresults(out.ROC.P,Pxnull,Pynull,[string '_P '],t);    
end


function [FA,pFA] = myAnova(X,Y)
         XG = cell(size(X));
         XG(find(X==-1)) = {'AA'};
         XG(find(X==1)) = {'BB'};
         XG(find(X==0)) = {'AB'};
         [~,T] = anova1(Y,XG,'off');
         FA = T{2,5};
         pFA  = T{2,6};
end

function M = myBootstrapRegressSingleGene(A,B,Model,type)
         Est = A(:,3)';
         for i=1:1:3
             Est = myregressSingleGeneLOO([A(:,1:2) Est'],B,Model,type);
         end
         [~,~,~,~,M,~,~,~] = plsregress([A(:,1:2) Est'],B,size(A,2));
         M = M(end,:);
end


% function [Est] = myBootstrapRegressSingleGeneLOO(A,B,Model,type)
%          Est = A(:,3)';
%          for i=1:1:3
%              Est = myregressSingleGeneLOO([A(:,1:2) Est'],B,Model,type);
%          end
% end

function [Est] = myregressSingleGeneLOO(A,B,Model,type)
         n = size(A,1);
         Est = zeros(1,n);
         Ind = (1:n);
         parfor i=1:n
            TrInd = setdiff(Ind,i);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            TestB = B(i,:);
            [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
            [~,~,Est(i),~] = getDistance(Model,[],TestB,type,TrM(end,:));
         end
end

function [out] = myROCtests(X,Y,msH,msM,msP,change)
         % Perform Homozygote ROC
         index = find(X);
         XH = X(index);
         XH(find(XH == -1)) = 0;
         YH = Y(index);
         [xH,yH,aucH,seH] = roc_calc(~logical(XH),YH*msH);
         if aucH < 0.5 && change
            msH = -1;
            [xH,yH,aucH,seH] = roc_calc(~logical(XH),YH*msH);
         end
         out.H.x = xH;
         out.H.y = yH;
         out.H.auc = aucH;
         out.H.se = seH;
         % Perform -1 ROC
         XM = X;
         XM(find(XM==0)) = -1;
         XM(find(XM == -1)) = 0;
         YM = Y;
         [xM,yM,aucM,seM] = roc_calc(~logical(XM),YM*msM);
         if aucM < 0.5 && change
            msM = -1;
            [xM,yM,aucM,seM] = roc_calc(~logical(XM),YM*msM);
         end
         out.M.x = xM;
         out.M.y = yM;
         out.M.auc = aucM;
         out.M.se = seM;
         % Perform +1 ROC
         XP = X;
         XP(find(XP==0)) = 1;
         XP(find(XP == -1)) = 0;
         YP = Y;
         [xP,yP,aucP,seP] = roc_calc(~logical(XP),YP*msP);
         if aucM < 0.5 && change
            msP = -1;
            [xP,yP,aucP,seP] = roc_calc(~logical(XP),YP*msP);
         end
         out.P.x = xP;
         out.P.y = yP;
         out.P.auc = aucP;
         out.P.se = seP;
end

function fig = plotROCresults(in,samplex,sampley,string,t)
         fig = figure;hold on;
         T = [string ' AUC:' num2str(in.auc) ' P:' num2str(in.pauc) ' Pperm:' num2str(in.ppermauc)];
         title(T);
         xedges = linspace(0,1,64); yedges = linspace(0,1,64);
         histmat = hist2(samplex(:), sampley(:), xedges, yedges);
         pcolor(xedges,yedges,histmat'/t); set(gca,'clim',[0 1]);colorbar ; axis square tight ;
         plot(in.x,in.y,'w-','Linewidth',2);
         xl = (0:0.01:1);
         yl = (0:0.01:1);
         plot(xl,yl,'w.','Linewidth',1.5);
end

