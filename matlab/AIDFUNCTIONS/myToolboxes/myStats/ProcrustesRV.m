function out = ProcrustesRV(X1,X2,t)
         if nargin < 3, t = 0; end
         % Rv coefficient between two equi-sized matrices
         X = [X1(:) X2(:)];
         num_RV = X'*X;
         la_diag=diag(num_RV);
         nd=length(la_diag);
         den_RV=(repmat(la_diag',nd,1).*repmat(la_diag,1,nd)).^(1/2);
         c = num_RV./den_RV;
         c = c(1,2);
         out.RV = c;
         if t==0, return; end
         n = numel(X1);
         CCount = false(t,1);
         tic;
         nc = size(X1,1);
         parfor i=1:t
                index = randperm(n);
                X1for = X1;
                X1for(1:n) = X1for(index);
                outfor = RV(X1for,X2,0);
                CCount(i) = outfor.RV>=c;
         end
         toc;
         out.permRV = sum(CCount)/t;         
end