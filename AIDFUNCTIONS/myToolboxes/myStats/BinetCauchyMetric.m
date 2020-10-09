function out = BinetCauchyMetric(A,n)
         if nargin < 2, n = length(A); end
         A = A(1:n);
         out = sqrt(1-prod(cos(A).^2));
end