function [LabelTotal,MASK] =HierarchicalFacialSegmentationv4(SimilarityMatrix,nLevels,type,nRuns,minPercLM)
    % This function creates sub-modules of the face through iterative runs of
    % the spectral clusteinrg (Ng et al.)
    % INPUT
    % SimilarityMatrix: n x n, with n the number of 3D landmarks in a face.
    % (this is the RV matrix as obtained from getRVmatrix)
    % n_levels: the number of levels in the hierarchical segmentation tree
    % Type: the type of the laplacian transformation applied to the
    % SimilarityMatrix, default is Weiss
    % runs: the amount of runs in kmeans++, to suppress random behaviour,
    % default = 50;
    % Output: LabelTotal, L by n matrix, giving cluster labels per level
    % (different levels in different rows) for each of the n data points (each
    % collumn)

    % check input 
    if nargin<3,type='weiss';end
    if nargin<4, nRuns = 50; end
    warning off;
    nCluster = 1;
    nLM = size(SimilarityMatrix,1);
    LabelTotal = ones(nLevels+1,nLM);
    h=waitbar(0,'please wait');
    for level=1:nLevels
        clustercount = 0;  
        for cluster=1:nCluster
            disp([datestr(now) ' Level ', num2str(level), ' start', ' cluster ', num2str(cluster), ' start']);
            pos = find(LabelTotal(level,:)==cluster);
            percLM = 100*(length(pos)/nLM);
            if percLM>minPercLM
               A = SimilarityMatrix(pos,pos);
               label = mySpectralClustering(A,type,nRuns);
               label(label==1)=clustercount+1;
               label(label==2)=clustercount+2;
               LabelTotal(level+1,pos)=label;
            else
               LabelTotal(level+1:end,pos) = nan;
            end
            clustercount = clustercount+2;
            disp([datestr(now) ' Level ', num2str(level), ' done', ' cluster ', num2str(cluster), ' done']); 
        end
        nCluster = clustercount;
        waitbar(level/nLevels,h);
    end
    close(h);
    HI = HierarchicalInterface;
    HI.nL = nLevels+1;
    MASK = ones(1,HI.nLC);
    for i=1:1:HI.nLC
        [l,c]= Ind2LC(HI,i);
        if isempty(find(LabelTotal(l,:)==c)), MASK(i)=0; end %#ok<EFIND>
    end
    warning on;
end

function [C_opt,cost]=mySpectralClustering(W,type,runs)
%   Executes spectral clustering algorithm
%    the default (weiss) was defined by Ng. et al in 'On Spectral Clustering: Analysis and an algorithm';
%   the spectral clustering is performed on the adjacency matrix W and returns the 2 cluster
%   indicator vectors as columns in C.
%   
%   'W' - Adjacency matrix, needs to be square
%   'k' - Number of clusters to look for
%    Spectral Clustering Algorithm Type: Normalized according to Jordan and Weiss (2002)
%
%   References:
%   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering",
%     Statistics and Computing 17 (4), 2007
    n = size(W,1);
    W(eye(n)==1) = 0;
    switch lower(type)
        case 'ratiocute'
            D = sum(W);
            L = diag(D)-W;
            [U,~] = eigs(L,2,'smallestabs');
            U = U(:,2);
            cost = 0;
            C_opt = double(U>0)+1;
            return;
        case 'ncute'
            D = sum(W);
            D = diag(D.^(-0.5));
            I = D*W*D;
            [U,~] = eigs(I,2);
            U = U(:,2);
            cost = 0;
            C_opt = double(U>0)+1;
            return;
        case 'symmetric laplacian'
            D = sum(W);
            D = diag(D.^(-0.5));
            I = D*W*D;
            [U, values] = eig(I,'vector');
            [~,ind] = sort(values,'descend');
            U = D*U(:,ind(2:3));
        case 'weiss' % DEFAULT
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
            k = 2;
            [U,~] = eigs(L, k, diff);
            % in case of the Jordan-Weiss algorithm, we need to normalize
            % the eigenvectors row-wise
            U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2))); 
        otherwise
    end
    cost = +Inf;
    k = 2;
    for i=1:runs
        [C,~,c] = myKmeans(U',k);
        %disp(num2str(c));
        if c < cost
            C_opt = C;
            cost = c;
        end
    end
end

function [L,C,cost] = myKmeans(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.
    L = [];
    L1 = 0;
    while length(unique(L)) ~= k

        % The k-means++ initialization.
        C = X(:,1+round(rand*(size(X,2)-1)));
        L = ones(1,size(X,2));
        for i = 2:k
            D = X-C(:,L);
            D = cumsum(sqrt(dot(D,D,1)));
            if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
            C(:,i) = X(:,find(rand < D/D(end),1));
            [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
        end

        % The k-means algorithm.
        while any(L ~= L1)
            L1 = L;
            for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
            [m,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
            cost = -1*sum(m);
        end
    end
end



