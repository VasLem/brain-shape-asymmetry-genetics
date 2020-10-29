function out = ProjectionMetric(A,n)
         if nargin < 2, n = length(A); end
         A = A(1:n);
         m = length(A);
         out = sqrt(m-sum(cos(A).^2));
end