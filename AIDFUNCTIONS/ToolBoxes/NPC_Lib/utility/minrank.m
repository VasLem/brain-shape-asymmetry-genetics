function [r] = minrank(x)
%minRANK Compute the ranks of a sample, adjusting for ties.
%   [R] = MINRANK(X) computes the ranks of the values in the
%   vector X.  If any X values are tied, MINRANK computes their MIN
%   rank. See also the original matlab function tiedrank
%
%livio finos

[n m]=size(x);
if nargin < 2
    flag = false;
end
[sx, rowidx] = sort(x);

for i=1:m
    ranks = (1:n)';
if isa(x(:,i),'single')
   ranks = single(ranks);
end

% Adjust for ties
tieloc = find(~diff(sx(:,i)));
    tieadj = 0;
while (length(tieloc) > 0)
    tiestart = tieloc(1);
    ntied = 1 + sum(sx(tiestart,i) == sx(tiestart+1:end,i));
    ranks(tiestart:tiestart+ntied-1) = tiestart;
    tieloc(1:ntied-1) = [];
end
r(rowidx(:,i),i) = ranks;
end
