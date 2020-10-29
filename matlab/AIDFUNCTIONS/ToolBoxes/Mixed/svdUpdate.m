function [U, s, V] = svdUpdate( U, s, V, C ) 
%
% Incremental SVD 
%  
% -------------------------------------------------------------------------
%  [U, s, V] = svdUpdate( U, s, V, C );
%
%
%   Input:  DOUBLE U (p by r)  collumn space
%           DOUBLE s (r by 1)  singular values
%           DOUBLE V (q by r)  row space
%           DOUBLE C (p by c)  additional collumns of data 
%
%  Output:  DOUBLE U (p by ?)  collumn space
%           DOUBLE s (? by 1)  singular values
%           DOUBLE V (q by ?)  row space 
%             
% -------------------------------------------------------------------------
%
% An implementation of the incremental SVD described in section 3 of
% Matthew Brands 2002 ECCV [1] paper. 
%  
%  [1] Brand, M.E., "Incremental Singular Value Decomposition of
%  Uncertain Data with Missing Values", European Conference on Computer
%  Vision (ECCV), Vol 2350, pps 707-720, May 2002 
%  
% An example:   
%
% $$$    %
% $$$    % generate a matrix M
% $$$    %
% $$$ 
% $$$    M = rand(50, 10);
% $$$ 
% $$$    %
% $$$    % compute the SVD of M and get the rank-r approximation
% $$$    %
% $$$ 
% $$$   [U S V] = svd( M, 0 );
% $$$ 
% $$$   r = min( size(U,1), size(V,2) );
% $$$ 
% $$$   U = U(:,1:r);
% $$$   s = diag(S(1:r,1:r)); 
% $$$   V = V(:,1:r);
% $$$ 
% $$$ 
% $$$   %
% $$$   % update the SVD 100 times
% $$$   %
% $$$ 
% $$$   for ii = 1:100
% $$$ 
% $$$ 	
% $$$     %
% $$$     % generate some random collumns to add to M  
% $$$     %
% $$$ 
% $$$     C = rand( size(U,1), 20); 
% $$$     
% $$$     %
% $$$     % compute the SVD of [M|C]
% $$$     %
% $$$      
% $$$     [U, s, V] = svdUpdate( U, s, V, C ); 
% $$$  
% $$$     %
% $$$     % update M
% $$$     %
% $$$ 
% $$$     M = [M C];
% $$$ 
% $$$     %
% $$$     % compute the reprojection error (should be around 1e-10 )
% $$$     %
% $$$ 	 
% $$$     norm(M - U*diag(s)*V') 
% $$$     
% $$$   end;
%
% -------------------------------------------------------------------------
%  
% Author:      Nathan Faggian
% E-mail:      nathanf@mail.csse.monash.edu.au
% URL:         http://www.csse.monash.edu.au/~nathanf
%    
   
  
%
% add a column C or matrix which is p*c
%

c = size(C, 2);
r = size(s, 1);
q = size(V, 1);

%
% compute the projection of C onto the orthogonal subspace U
%

L = U'*C;

%
% compute the component of C orthogonal that is orthogonal to the subspace U
%

H = C - U*L;

%
% compute an orthogonal basis of H and the projection of C onto the
% subspace orthogonal to U
%

[J, K] = qr(H,0);
 
%
% compute the center matrix Q
%
 
Q = [      diag(s),  L;
        zeros(c,r),  K   ];

%
% compute the SVD of Q 
%

[Uu, Su, Vu] = svd(Q, 0);

%
% compute the updated SVD of [M,C]
%

orth = U(:,1)'*U(:,end);

U = [U, J] * Uu;

s = diag(Su);

V = [   V          , zeros(q, c); ...
        zeros(c, r), eye(c)     ] * Vu;

%
% compact the new SVD
%

r = min( size(U,1), size(V,2) );

U = U(:,1:r);
s = s(1:r);
V = V(:,1:r);






