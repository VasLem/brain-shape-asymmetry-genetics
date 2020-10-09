function out = multiplePlsRegress(A,B,t,type)         
         n = size(A,1);
         nA = size(A,2);
         ncomp = min(n,nA);
         if nargin <4, type = 'value'; end
         if nargin < 3, t = 0; end
         % multiple regression to obtain multiple R-Squared
         [~,~,~,~,M,pctvar,~,stats] = plsregress(A,B,ncomp);
         R2 = sum(pctvar(2,:));
         R2S = shrinkR2(R2,n,nA);
         [LocalS,LocalE] = getLocalStatistic(M,stats.Yresiduals,B,type);
         out.R2 = R2;
         out.R2S = R2S;
         out.LocalE = LocalE;
         out.LocalR2 = LocalS;
         if t==0, return; end
         LocalSCount = false(length(LocalS),t);
         R2Count = false(t,1);
         %tic;
         parfor i=1:t % permutation loop
                ind = randperm(n);
                Bfor = B(ind,:);  %#ok<PFBNS>
                [~,~,~,~,Mfor,pctvarfor,~,statsfor] = plsregress(A,Bfor,ncomp);
                R2for = sum(pctvarfor(2,:));
                R2Count(i) = (R2for>=R2);
                [LocalSfor] = getLocalStatistic(Mfor,statsfor.Yresiduals,Bfor,type);
                LocalSCount(:,i) = (LocalSfor >= LocalS);
         end % end permutation loop
         %toc;
         out.pR2 = sum(R2Count)/t;
         out.pLocalR2 = sum(LocalSCount,2)/t;
end

function [LocalS,LocalE] = getLocalStatistic(M,R,B,type)
         [n,nB] = size(B);
         nA = size(M,1)-1;
         P = B-R;
         A = repmat(mean(B),n,1);
         SST = sum((B-A).^2);
         SSR = sum((P-A).^2);
         switch lower(type)
             case {'value' 'pc'}
                 LocalS = SSR./SST;
                 if nargout==1, return; end
                 LocalE = M(2:end,:);
             case 'shape'
                 SSR = reshape(SSR,3,nB/3);
                 SSR = sum(SSR);
                 SST = reshape(SST,3,nB/3);
                 SST = sum(SST);
                 LocalS = SSR./SST;
                 if nargout==1, return; end
                 LocalE = zeros(4,nB/3,nA);
                 for i=1:1:nA
                     LocalE(1:3,:,i) = reshape(M(i+1,:),3,nB/3);
                     LocalE(4,:,i) = sqrt(sum(LocalE(1:3,:,i).^2));
                 end
             otherwise
                 error('unknown type');
         end        
end


function out = shrinkR2(R2,n,nA)
         out = 1-(1-R2)*((n-1)/(n-nA-1));
end

