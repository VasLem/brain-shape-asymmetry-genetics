function [out] = myPlsRegress(A,B,t,p,type)         
         n = size(A,1);
         nA = size(A,2);
         nB = size(B,2);
         ncomp = min(n,nA);
         if nargin <5, type = 'value'; end
         if nargin < 4, p = (1:nA); end
         if nargin < 3, t = 0; end
         [~,~,~,~,M,pctvar,~,stats] = plsregress(A,B,ncomp);
         out.M = M;
         out.pctvar = pctvar;
         out.Tpctvar = sum(out.pctvar,2);
         out.Stats = stats;
         [out.LocalE,out.LocalS,out.GlobalS] = getRegressionStatistic(M,nA,nB,type);
         if t==0, return; end
         f = statusbar('Permuting');
         drawnow;
         for k=1:1:length(p);
             % creating reduced model
             index = setdiff((1:nA),p(k));
             AR = A(:,index);
             [~,~,~,~,~,~,~,statsR] = plsregress(AR,B,ncomp-1);
             E = statsR.Yresiduals;
             Afor = A(:,p(k));
             % getting reference test statistics
             [~,~,~,~,M,pctvar,~,~] = plsregress(Afor,E,1);
             R2 = pctvar(2); 
             [LocalE,LocalS,GlobalS] = getRegressionStatistic(M,1,nB,type);
             LocalSCount = false(length(LocalS),t);
             R2Count = false(t,1);
             GlobalSCount = false(t,1);
             parfor i=1:t % permutation loop
                    ind = randperm(n);
                    Bfor = E(ind,:);  %#ok<PFBNS>
                    [~,~,~,~,Mfor,pctvarfor,~,~] = plsregress(Afor,Bfor,1);
                    R2for = pctvarfor(2);
                    R2Count(i) = (R2for>=R2);
                    [~,LocalSfor,GlobalSfor] = getRegressionStatistic(Mfor,1,nB,type);
                    LocalSCount(:,i) = (LocalSfor >= LocalS);
                    GlobalSCount(i) = GlobalSfor >= GlobalS;
              end % end permutation loop
              PermTest(p(k)).R2 = R2;
              PermTest(p(k)).pR2 = sum(R2Count)/t;
              PermTest(p(k)).LocalE = LocalE;
              PermTest(p(k)).LocalS = LocalS;
              PermTest(p(k)).pLocalS = sum(LocalSCount,2)/t;
              PermTest(p(k)).GlobalS = GlobalS;
              PermTest(p(k)).pGlobalS = sum(GlobalSCount)/t;
              statusbar(k/length(p),f);
              drawnow;
         end 
         delete(f);
         out.PermTest = PermTest;
end


