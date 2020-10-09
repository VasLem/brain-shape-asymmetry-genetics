function out = getModelMemory(obj,G,rank,startN,nr)
% out = getModelMemory(obj,G,rank,startN,n)
% Function to obtain a model in a more memory friendly way, takes more
% time.
% INPUT
% obj = PCA object
% G = data structure or matrix
% rank = maximum number of pc's to retain (please limit in order to use the
% memory efficiency)
% startN = the number of start measurements to use
% n = the number of measurements to use in each update.
% OUTPUT
% out = PCA object with eigenvectors and eigenvalues
%
%
% Created by Peter Claes
    if nargout == 1,obj = clone(obj);out = obj;end
    if nargin < 4,startN = 10;end
    if nargin < 5, nr = 1; end
    if nargin < 3, rank = 200; end
    G = getData(obj,G);
    if isempty(obj.AvgVec), getAverage(obj,G); end
    n = size(G,2);
    %nanIndex = find(isnan(D));
    % Initial Model;
    [G,D] = selectData(G,startN);
    if obj.Centering, 
       [U,S,V] = svd(double(D)-repmat(obj.AvgVec,1,size(D,2)),'econ');
    else
        [U,S,V] = svd(double(D),'econ');
    end
    r = min([size(U,1) rank size(V,2)]);
    U = U(:,1:r);
    S = diag(S(1:r,1:r)); 
    V = V(:,1:r);
    f = statusbar('Memory efficient Model Building');
    while ~isempty(G)
        [G,D] = selectData(G,nr);
        if obj.Centering
            [U, S, V] = svdUpdate( U, S, V, double(D)-repmat(obj.AvgVec,1,size(D,2)));
        else
            [U,S,V] = svd(double(D),'econ');
        end
        if size(U,2)>rank
           U = U(:,1:rank);
           S = S(1:rank);
           V = V(:,1:rank);
        end
        clear D;
        statusbar((1-size(G,2)/n),f);
    end
    delete(f);
    % EigVec
    obj.EigVec = U;
    % Coefficients
    obj.Tcoeff = repmat(S',n,1).*V;
    % EigenValues
    obj.EigVal = (S./sqrt(n-1)).^2;
    %   
    if n<rank
       obj.EigVal(n) = [];
       obj.EigVec(:,n) = [];
       obj.Coeff(:,n) = [];
    end

end

function [GF,GP] = selectData(G,nr)
         try
             GP = G(:,1:nr);
             GF = G(:,nr+1:end);
         catch %#ok<CTCH>
             GP = G;
             GF = [];
         end
end

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

    % add a column C or matrix which is p*c
    c = size(C, 2);
    r = size(s, 1);
    q = size(V, 1);
    % compute the projection of C onto the orthogonal subspace U
    L = U'*C;
    % compute the component of C orthogonal that is orthogonal to the subspace U
    H = C - U*L;
    % compute an orthogonal basis of H and the projection of C onto the
    % subspace orthogonal to U
    [J, K] = qr(H,0);
    % compute the center matrix Q
    Q = [      diag(s),  L;
            zeros(c,r),  K   ];
    % compute the SVD of Q 
    [Uu, Su, Vu] = svd(Q, 0);
    % compute the updated SVD of [M,C]
    orth = U(:,1)'*U(:,end);
    U = [U, J] * Uu;
    s = diag(Su);
    V = [   V          , zeros(q, c); ...
            zeros(c, r), eye(c)     ] * Vu;
    % compact the new SVD
    r = min( size(U,1), size(V,2) );
    U = U(:,1:r);
    s = s(1:r);
    V = V(:,1:r);
    
end