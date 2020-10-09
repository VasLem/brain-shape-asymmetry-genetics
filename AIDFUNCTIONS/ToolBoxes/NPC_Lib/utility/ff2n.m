function x = ff2n(n)
%ff2n   Two-level full-factorial design.
%   X = ff2n(N) creates a two-level full-factorial design, X.
%   N is the number of columns of X. The number of rows is 2^N.

%   B.A. Jones 2-17-95
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.8 $  $Date: 2002/01/17 21:30:32 $

rows = 2.^(n);
ncycles = rows;
x = zeros(rows,n);

for k = 1:n
   settings = (0:1);
   ncycles = ncycles/2;
   nreps = rows./(2*ncycles);
   settings = settings(ones(1,nreps),:);
   settings = settings(:);
   settings = settings(:,ones(1,ncycles));
   x(:,n-k+1) = settings(:);
end
