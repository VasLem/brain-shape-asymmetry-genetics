function [C_opt,L,U]=spectral_clustering_weiss_ng(W, k)
%   spectral_clustering_weiss_ng Executes spectral clustering algorithm
%   defined by Ng. et al in 'On Spectral Clustering: Analysis and an algorithm';
%   the spectral clustering is performed on the adjacency matrix W and returns the k cluster
%   indicator vectors as columns in C.
%   If L and U are also called, the (normalized) Laplacian and
%   eigenvectors will also be returned.
%
%   'W' - Adjacency matrix, needs to be square
%   'k' - Number of clusters to look for
%    Spectral Clustering Algorithm Type: Normalized according to Jordan and Weiss (2002)
%
%   References:
%   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering",
%     Statistics and Computing 17 (4), 2007


% calculate degree matrix
degs = sum(W, 2); % sum of the elements along the row
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);

% compute unnormalized Laplacian
L = D - W;


% compute normalized Laplacian 
% avoid dividing by zero
degs(degs == 0) = eps;
% calculate D^(-1/2)
D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));

% calculate normalized Laplacian
L = D * L * D;


% compute the eigenvectors corresponding to the k smallest
% eigenvalues
diff = eps;
[U,~] = eigs(L, k, diff);

% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));

% now use the k-means algorithm to cluster U row-wise
% C will be a n-by-1 matrix containing the cluster number for
% each data point

% C = kmeans(U, k, 'start', 'cluster', ...
%     'EmptyAction', 'singleton'); % option 'plus' cannot run in Matlab
%     R2012b

%c_opt = -Inf;
c_opt = +Inf;
for i=1:50
    [C,~,c] = kmeans_Laurent_Sorber(U',k);
    disp(num2str(c));
    if c < c_opt
        C_opt = C;
        c_opt = c;
    end
end
% disp(num2str(c_opt))
end