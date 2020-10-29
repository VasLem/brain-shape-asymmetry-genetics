function [d,A]  = getSubspaceDistance(QF,QG,type)
         if nargin < 3, type = 'Projection Frobenius'; end
         %if nargin < 4, nA = +inf; end
         k1 = size(QF,2);
         k2 = size(QG,2);
         k = min(k1,k2);
         A = mySubspaceAngles(QF,QG);
         %if length(A)>nA, A = A(1:nA); k = nA; k1 = nA; k2 = nA; end
         switch lower(type)
             case 'projection frobenius'
                 d = sqrt(k1+k2-2*sum(cos(A).^2));
             case 'projection'
                 d = sqrt(k-sum(cos(A).^2));
             case 'procrustes'
                 d = 2*sqrt(sum(sin(A/2).^2));
             case 'mitteroecker'
                 d = sqrt(sum(log(cos(A)).^2));
             case 'average'
                 d = k/sum(cos(A));
         end
end