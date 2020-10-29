function out = SubspaceDStatistic(S1,S2,t)
         S1 = orth(S1);
         S2 = orth(S2);
         [DStat, A] = getSubspaceDistances(S1,S2,'projection');
         out.S1 = S1;
         out.S2 = S2;
         out.DStat = DStat;
         out.A = A;
         if t<=0, return; end
         StatCount = false(length(DStat),t);
         DStatpermC = zeros(t,length(DStat));
         disp('Random Subspaces');
         tic;
         parfor i=1:t
             R = randn(size(S1));
             R = orth(R);
             DStatR = getSubspaceDistances(S1,R,'projection');
             StatCount(:,i) = (DStatR<=DStat)';
             DStatpermC(i,:) = DStatR;
         end       
         toc;
         figure;
         set(gca,'Xlim',[1 length(DStat)]);
         hold on;
         for i=1:1:t
             plot(DStatpermC(i,:),'r-.');
             drawnow;   
         end
         plot(DStat,'b','linewidth',3);
         out.p = (sum(StatCount,2)/t)';
         out.DStatpermC = DStatpermC;
end