function out = ProcrustesMetric(A,n)
         if nargin < 2, n = length(A); end
         A = A(1:n);
         out = 2*sqrt(sum(sin(A/2).^2));
end