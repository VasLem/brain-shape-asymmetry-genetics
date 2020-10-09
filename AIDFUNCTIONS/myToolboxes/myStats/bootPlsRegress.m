function boot = bootPlsRegress(A,B,t,p)         
         if nargin < 4, p = (1:size(A,2)); end
         if nargin < 3, t = 0; end
         if t< 1000, t = 0; end
         n = size(A,1);
         nA = size(A,2);
         ncomp = min(n,nA);
         Tboot = zeros(1,t);
         parfor i=1:t
            % bootstrap
            ind = randsample((1:n),n,true);
            Afor = A;
            Bfor = B;
            Afor(:,p) = Afor(ind,p);
            Bfor(:,p) = Bfor(ind,p);
            [~,~,~,~,~,pctvarfor] = plsregress(Afor,Bfor,ncomp);
            Tboot(i) = sum(pctvarfor(2,:));
         end
         Tboot = sort(Tboot);
         boot.Tboot = Tboot;
         boot.upperci = Tboot(0.975*t);
         boot.lowerci= Tboot(0.025*t+1);
         boot.Avg = mean(Tboot);
end

% if nargin < 3, t = 1000; end
%          if t< 1000, t = 1000; end
%          n = size(A,1);
%          nA = size(A,2);
%          ncomp = min(n,nA);
%          [~,~,~,~,M,pctvar] = plsregress(A,B,ncomp);      
%          BR = B-[ones(n,1) A]*M;
%          denom = det(cov(B));
%          L = 1-det(cov(BR))/denom;
%          pctvarcount = false(2,ncomp,t);
%          Ls = zeros(1,t);
%          tic
%          parfor i=1:t
%              ind=randperm(n);
%              Afor = A(ind,:);
%              [~,~,~,~,Mfor,pctvarfor] = plsregress(Afor,B,ncomp);
%              pctvarcount(:,:,i) = pctvarfor>=pctvar;
%              BRfor = B-[ones(n,1) Afor]*Mfor;
%              Ls(i) = 1-det(cov(BRfor))/denom;
%          end
%          toc
%          pctvarP = (sum(pctvarcount,3)/t);
%          [sortLs]=sort(Ls);
%          upperci=sortLs(0.975*t);
%          Lsign = L>upperci;       