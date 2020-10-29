function out = MitteroeckerMetric(A,n)
         if nargin < 2, n = length(A); end
         A = A(1:n);
         out = sqrt(sum(log(cos(A)).^2));
end