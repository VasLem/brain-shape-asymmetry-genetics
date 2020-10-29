function PartialTest = partialPlsRegress(A,B,t,p,type)
         [A,B] = eliminateNAN(A,B);
         n = size(A,1);
         nA = size(A,2);
         nB = size(B,2);
         ncomp = min(n,nA);
         if nargin <5, type = 'value'; end
         if nargin < 4, p = (1:nA); end
         if nargin < 3, t = 1000; end
         %if t==0, PartialTest = []; return; end
         %f = statusbar('Partial Regr.');
       
         % multiple regression to obtain multiple R-Squared
         [~,~,~,~,~,pctvar] = plsregress(A,B,ncomp);
         R2Tot = sum(pctvar(2,:)); 
         drawnow;
         for k=1:1:length(p);
             disp(num2str(p(k)));
             % creating reduced models
             index = setdiff((1:nA),p(k));
             AR = A(:,index);
             [~,~,~,~,~,pctvar,~,statsR] = plsregress(AR,B,ncomp-1);
             E = statsR.Yresiduals;
             R2Red = sum(pctvar(2,:));
             R2Add = R2Tot-R2Red;
             Afor = A(:,p(k));
             [~,~,~,~,~,~,~,statsR] = plsregress(AR,Afor);
             Afor = statsR.Yresiduals;
             % getting reference test statistics
             [~,~,~,~,M,pctvar,~,stats] = plsregress(Afor,E,1);
             R2Part = pctvar(2);
             [LocalE,LocalS] = getLocalStatistic(M,stats.Yresiduals,E,type);
             PartialTest(p(k)).R2Tot = R2Tot; %#ok<*AGROW>
             PartialTest(p(k)).R2Red = R2Red; 
             PartialTest(p(k)).R2Add = R2Add;
             PartialTest(p(k)).R2Part = R2Part;
             PartialTest(p(k)).LocalE = LocalE;
             PartialTest(p(k)).LocalR2 = LocalS;
             PartialTest(p(k)).M = M;
             if t==0, continue; end
             LocalSCount = false(length(LocalS),t);
             R2PartCount = false(t,1);
              disp(['PLSR VAR' num2str(p(k))]);
             [path,ID] = setupParForProgress(t);
             parfor i=1:t % permutation loop
                    ind = randperm(n);
                    Bfor = E(ind,:);  %#ok<PFBNS>
                    [~,~,~,~,Mfor,pctvarfor,~,statsfor] = plsregress(Afor,Bfor,1);
                    R2for = pctvarfor(2);
                    R2PartCount(i) = (R2for>=R2Part);
                    [~,LocalSfor] = getLocalStatistic(Mfor,statsfor.Yresiduals,Bfor,type);
                    LocalSCount(:,i) = (LocalSfor >= LocalS);
                        parfor_progress;
                    %GlobalSCount(i) = GlobalSfor >= GlobalS;
             end % end permutation loop
             closeParForProgress(path,ID);
             PartialTest(p(k)).pR2Part = sum(R2PartCount)/t;
             PartialTest(p(k)).pLocalR2 = sum(LocalSCount,2)/t;
%              statusbar(k/length(p),f);
%              drawnow;
          
         end 
         %delete(f);
end


function [LocalE,LocalS] = getLocalStatistic(M,R,B,type)
         [n,nB] = size(B); 
         P = B-R;
         A = repmat(mean(B),n,1);
         SST = sum((B-A).^2);
         SSR = sum((P-A).^2);
         switch lower(type)
             case {'value' 'pc'}
                 LocalE = M(2,:);
                 LocalS = SSR./SST;
              case 'shape'
                 LocalE = zeros(4,nB/3);
                 LocalE(1:3,:) = reshape(M(2,:),3,nB/3);
                 LocalE(4,:) = sqrt(sum(LocalE(1:3,:).^2));
                 SSR = reshape(SSR,3,nB/3);
                 SSR = sum(SSR);
                 SST = reshape(SST,3,nB/3);
                 SST = sum(SST);
                 LocalS = SSR./SST;
             otherwise
                 error('unknown type');
         end        

end

