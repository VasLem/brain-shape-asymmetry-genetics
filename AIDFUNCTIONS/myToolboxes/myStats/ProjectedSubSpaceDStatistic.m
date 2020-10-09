function out = ProjectedSubSpaceDStatistic(XA,XB,Y,t,ncmp)
         [~,~,~,~,MA] = plsregress(XA,Y,ncmp);
         [~,~,~,~,MB] = plsregress(XB,Y,ncmp);
         SA = orth(MA(2:end,:)');
         SB = orth(MB(2:end,:)');
         [DStat, A, U, V] = getSubspaceDistances(SA,SB,'projection');
         out.SA = SA;
         out.SB = SB;
         out.DStat = DStat;
         out.A = A;
         out.U = U;
         out.V = V;
         if t<=0, return; end
         StatCount = false(length(DStat),t);
         DStatpermC = zeros(t,length(DStat));
         disp('permuting');
         tic;
         nT = size(XA,1);
         parfor i=1:t
             ind = randperm(nT);
             XAFor = XA(ind,:); %#ok<PFBNS>
             [~,~,~,~,MAFor] = plsregress(XAFor,Y,ncmp);
             SAFor = orth(MAFor(2:end,:)');
             DStatFor = getSubspaceDistances(SB,SAFor,'projection');
             StatCount(:,i) = (DStatFor<=DStat)';
             DStatpermC(i,:) = DStatFor;
         end
         toc;
         try
             figure;
             set(gca,'Xlim',[1 length(DStat)]);
             hold on;
             for i=1:1:t
                 plot(DStatpermC(i,:),'r-.');
                 drawnow;   
             end
             plot(DStat,'b','linewidth',3);
         catch
         end
         out.p = (sum(StatCount,2)/t)';
         out.DStatpermC = DStatpermC;
end