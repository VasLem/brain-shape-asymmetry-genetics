function [out,fig] = permSGNullAnalysis(A,B,Model,t,string,type)
         if nargin < 4, t = 100; end
         ms = 1;
         n = size(A,1);
         Est = zeros(1,n);
         Ind = (1:n);
         f = statusbar('LOO');drawnow;
         %type = 'Euclidean';
         for i=1:1:n
            TrInd = setdiff(Ind,i);
            TrA = A(TrInd,:);
            TrB = B(TrInd,:);
            TestB = B(i,:);
            [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
            [~,~,Est(i),~] = getDistance(Model,[],TestB,type,TrM(end,:));
            statusbar(i/n,f);drawnow;
         end
         X = A(:,end);Y = ms*Est;
         X(find(X == -1)) = 0;
         [x,y,auc,se] = roc_calc(~logical(X),Y);
         if auc <0.5
             ms = -1;
             n = size(A,1);
             Est = zeros(1,n);
             Ind = (1:n);
             f = statusbar('LOO');drawnow;
             %type = 'Euclidean';
             for i=1:1:n
                TrInd = setdiff(Ind,i);
                TrA = A(TrInd,:);
                TrB = B(TrInd,:);
                TestB = B(i,:);
                [~,~,~,~,TrM,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
                [~,~,Est(i),~] = getDistance(Model,[],TestB,type,TrM(end,:));
                statusbar(i/n,f);drawnow;
             end
             X = A(:,end);Y = ms*Est;
             X(find(X == -1)) = 0;
             [x,y,auc,se] = roc_calc(~logical(X),Y);  
         end
         out.GeneDP = Est;
         out.auc = auc;
         out.se = se;
         out.x = x;
         out.y = y;
         delete(f);
         disp('Permuting');
         n = size(A,1);
         sampleauc = zeros(1,t);
         samplex = zeros(n,t);
         sampley = zeros(n,t);
         parfor b=1:t
            ind = randperm(n);
            Afor = A;
            Afor(:,end) = Afor(ind,end);
            Bfor = B;       
            Estfor = zeros(1,n);
            for i=1:1:n
                TrInd = setdiff(Ind,i);
                TrA = Afor(TrInd,:);
                TrB = Bfor(TrInd,:);
                TestB = Bfor(i,:);
                [~,~,~,~,M,~,~,~] = plsregress(TrA,TrB,size(TrA,2));
                [~,~,Estfor(i),~] = getDistance(Model,[],TestB,type,M(end,:));
            end
            X = Afor(:,end);Y = ms*Estfor;
            X(find(X == -1)) = 0;
            [samplex(:,b),sampley(:,b),sampleauc(b)] = roc_calc(~logical(X),Y);
         end
         out.ppermauc = sum(sampleauc>=auc)/t;
         out.pauc = normpdf((auc-0.5)/se,0,1);         
%          figure;hold on;
%          plot(x,y,'b-','LineWidth',2);
%          for i=1:10:t
%              plot(samplex(:,i),sampley(:,i),'r-');
%          end 
%          title(string);
         fig = figure;hold on;
         T = [string ' AUC:' num2str(out.auc) ' P:' num2str(out.pauc) ' Pperm:' num2str(out.ppermauc)];
         title(T);
         xedges = linspace(0,1,64); yedges = linspace(0,1,64);
         histmat = hist2(samplex(:), sampley(:), xedges, yedges);
         pcolor(xedges,yedges,histmat'/t); set(gca,'clim',[0 1]);colorbar ; axis square tight ;
         plot(out.x,out.y,'w-','Linewidth',2);
         xl = (0:0.01:1);
         yl = (0:0.01:1);
         plot(xl,yl,'w.','Linewidth',1.5);
end


%          sampleauc = sort(sampleauc);
%          out.upper = sampleauc(0.975*t);
%          out.lower= sampleauc(0.025*t+1);
%          out.avg = mean(sampleauc);
%          out.sampleAUC = sampleauc;
%          samplex = sort(samplex,2);
%          out.samplex = samplex;
%          out.avgx = mean(samplex,2);
%          out.upperx = samplex(:,0.975*t);
%          out.lowerx = samplex(:,0.025*t+1);
%          sampley = sort(sampley,2);
%          out.sampley = sampley;
%          out.avgy = mean(sampley,2);
%          out.uppery = sampley(:,0.975*t);
%          out.lowery = sampley(:,0.025*t+1);


% if nargin < 4, t = 100; end
%          n = size(A,1);
%          sampleauc = zeros(1,t);
%          samplex = zeros(n,t);
%          sampley = zeros(n,t);
%          Ind = (1:n);
%          parfor b=1:t
%             ind = randperm(n);
%             Afor = A;
%             Afor(:,end) = Afor(ind,end);
%             Bfor = B;
%             [~,~,~,~,M,~,~,~] = plsregress(Afor,Bfor,size(Afor,2));
%             Estfor = zeros(1,n);
%             for i=1:1:n
%                 TestB = Bfor(i,:);
%                 [~,~,Estfor(i),~] = getDistance(Model,[],TestB,'Mahalanobis',M(end,:));
%             end
%             X = Afor(:,end);Y = Estfor;
%             X(find(X == -1)) = 0;
%             [samplex(:,b),sampley(:,b),sampleauc(b)] = roc_calc(~logical(X),Y);
%          end
%          sampleauc = sort(sampleauc);
%          out.upper = sampleauc(0.975*t);
%          out.lower= sampleauc(0.025*t+1);
%          out.avg = mean(sampleauc);
%          out.sampleAUC = sampleauc;
%          samplex = sort(samplex,2);
%          out.samplex = samplex;
%          out.avgx = mean(samplex,2);
%          out.upperx = samplex(:,0.975*t);
%          out.lowerx = samplex(:,0.025*t+1);
%          sampley = sort(sampley,2);
%          out.sampley = sampley;
%          out.avgy = mean(sampley,2);
%          out.uppery = sampley(:,0.975*t);
%          out.lowery = sampley(:,0.025*t+1);
%          figure;hold on;
%          plot(out.avgx,out.avgy,'r-');
%          plot(out.upperx,out.uppery,'r.');
%          plot(out.lowerx,out.lowery,'r.');     