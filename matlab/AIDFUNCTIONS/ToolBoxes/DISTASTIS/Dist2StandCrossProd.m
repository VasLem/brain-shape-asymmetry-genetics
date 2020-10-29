function [rS,S] = Dist2StandCrossProd(D,m)
% Conversion of distance matrices to standarised Cross product matrices
     nfaces = size(D,1);
     njuges = size(D,3);
     if nargin < 2,m=ones(nfaces,1)./nfaces;end% facial weigth, equally important here
     cent=eye(nfaces,nfaces)-ones(nfaces,1)*m';
     Dsum=zeros(nfaces,nfaces);
     rS = D;
     S = D;
     f = waitbar(0,'CONVERTING');
     for k=1:njuges;
         Dk = D(:,:,k);
         Dsum = Dsum+Dk;
         % transform D -> SCP 
         rSk=-(1/2)*cent*Dk*cent;
         %[PS,lS]=eigen((rSk+rSk')/2);
         [~,lS] = svds(double((rSk+rSk')/2));% FASTER
         % normaliz by first eig
         Sk=rSk.*(lS(1).^(-1));% normalize
         rS(:,:,k) = rSk;
         S(:,:,k) = Sk;
         waitbar(k/njuges,f);
     end
     close(f);
end


% tic;[u,s]=svds(test,1);toc
% tic;[PS,lS]=eigen(test);toc;