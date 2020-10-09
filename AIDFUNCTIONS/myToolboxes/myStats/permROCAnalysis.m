function [out,fig] = permROCAnalysis(X,Y,t,string)
         if nargin < 3, t = 100; end
         ms = 1;
         [x,y,auc,se] = roc_calc(X,ms*Y);
         if auc <0.5
             ms = -1;
             [x,y,auc,se] = roc_calc(X,ms*Y);  
         end
         out.auc = auc;
         out.se = se;
         out.x = x;
         out.y = y;
         disp('Permuting');
         n = size(X,1);
         sampleauc = zeros(1,t);
         samplex = zeros(n,t);
         sampley = zeros(n,t);
         parfor b=1:t
            ind = randperm(n);   
            Yfor = Y(ind); %#ok<*PFBNS>
            [samplex(:,b),sampley(:,b),sampleauc(b)] = roc_calc(X, ms*Yfor);
         end
         out.ppermauc = sum(sampleauc>=auc)/t;
         out.pauc = normpdf((auc-0.5)/se,0,1);         
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