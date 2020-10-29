function CL = mySpecClustering(A,k)
    %perform Spectral Clustering
    %calculate degree matrix
    degs = sum(A, 2);
    D    = sparse(1:size(A, 1), 1:size(A, 2), degs);
    % compute unnormalized Laplacian
    L = D - A;
    % compute normalized Laplacian
    % avoid dividing by zero
    degs(degs == 0) = eps;
    % calculate inverse of D
    D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
    % calculate normalized Laplacian
    L = D * L;
    % Compute all the k eigenvector
    diff   = eps;
    [U, E] = eigs(L, k, diff);
    % k means clustering
    CL = kmeans(U, k, 'start', 'cluster','EmptyAction', 'singleton');
end