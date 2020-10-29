function [d, A, U, V] = getSubspaceDistances(QF,QG,type)
         if nargin < 3, type = 'Projection Frobenius'; end
         k1 = size(QF,2);
         k2 = size(QG,2);
         k = min(k1,k2);
         [A, U, V] = mySubspaceAngles(QF,QG);
         %if length(A)>nA, A = A(1:nA); k = nA; k1 = nA; k2 = nA; end
         nA = length(A);
         d = zeros(1,nA);
         for i=1:1:nA
             %k = i; k1 = i; k2 = i;
             k= i;
             switch lower(type)
                 case 'projection frobenius'
                     d(i) = sqrt(k1+k2-2*sum(cos(A(1:i)).^2));
                 case 'projection'
                     d(i) = sqrt(k-sum(cos(A(1:i)).^2));
                 case 'procrustes'
                     d(i) = 2*sqrt(sum(sin(A(1:i)/2).^2));
                 case 'mitteroecker'
                     d(i) = sqrt(sum(log(cos(A(1:i))).^2));
                 case 'average'
                     d(i) = k/sum(cos(A(1:i)));
                 case 'krzanowski'
                     d(i) = k-sum(cos(A(1:i)).^2);
             end
         end
end