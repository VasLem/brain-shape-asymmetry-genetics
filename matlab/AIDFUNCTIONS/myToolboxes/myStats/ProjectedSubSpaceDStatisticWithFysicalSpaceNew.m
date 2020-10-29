function out = ProjectedSubSpaceDStatisticWithFysicalSpaceNew(XA,Y,t,ncmp)
         [~,~,~,~,MA] = plsregress(XA,Y,size(XA,2));
         SA = orth(MA(2:end,:)');
         %SB = orth(MB(2:end,:)');
         SB = eye(size(Y,2),ncmp);
         [DStat, A, U, V] = getSubspaceDistances(SA,SB);
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
             [~,~,~,~,MAFor] = plsregress(XAFor,Y,size(XAFor,2));
             SAFor = orth(MAFor(2:end,:)');
             DStatFor = getSubspaceDistances(SB,SAFor);
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

function [D, A, U, V] = getSubspaceDistances(SA,SB)
         [Ys,s,Zs] = svd(SA'*SB,0);
         % principal vectors
         U = SA*Ys;
         V = SB*Zs;
         % principal angles
         A = acos(diag(s));
         nA = length(A);
         % subspace distances
         D = zeros(1,nA);
         for i=1:1:nA
             k = i;
             D(i) = sqrt(k-sum(cos(A(1:i)).^2));
         end
end