function out = DistastisSpaceCorrelation(X1,X2,type1,type2,t)
         if nargin < 3, t = 0; end
         if isempty(type2), type2 = type1;end
         [n,~] = size(X1);
         D(:,:,1) = squareform(pdist(X1,type1));
         D(:,:,2) = squareform(pdist(X2,type2));
         % Conversion of distance matrices to standarised Cross product matrices
           m = ones(n,1)./n;
           cent = eye(n,n)-ones(n,1)*m';
           rS = zeros(size(D));
           S = zeros(size(D));
           for k=1:1:2
               rSk=-(1/2)*cent*D(:,:,k)*cent;
               [~,lS]=eigen(rSk);
               % normaliz by first eig
               Sk=rSk.*(lS(1).^(-1));% normalize
               rS(:,:,k) = rSk;
               S(:,:,k) = Sk;
           end
         % Getting Rv coefficient
           out = RV(S(:,:,1),S(:,:,2),t);
end


% function out = DistastisSpaceCorrelation(X1,X2,type,t)
%          if nargin < 3, t = 0; end
%          [n,~] = size(X1);
%          if strcmpi(type,'Mahalanobis');
%            disp('Mahalanobis');
%          % Creation of Distance Matrices
%            D(:,:,1) = squareform(pdist(X1,'seuclidean'));
%            D(:,:,2) = squareform(pdist(X2,'seuclidean'));
%          else
%            disp('Euclidean');  
%            D(:,:,1) = squareform(pdist(X1));
%            D(:,:,2) = squareform(pdist(X2));  
%          end  
%          % Conversion of distance matrices to standarised Cross product matrices
%            m = ones(n,1)./n;
%            cent = eye(n,n)-ones(n,1)*m';
%            for k=1:1:2
%                rSk=-(1/2)*cent*D(:,:,k)*cent;
%                [~,lS]=eigen(rSk);
%                % normaliz by first eig
%                Sk=rSk.*(lS(1).^(-1));% normalize
%                rS(:,:,k) = rSk;
%                S(:,:,k) = Sk;
%            end
%          % Getting Rv coefficient
%            out = RV(S(:,:,1),S(:,:,2),t);
% end