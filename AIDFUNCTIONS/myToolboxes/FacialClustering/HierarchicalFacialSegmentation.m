function LabelTotal=HierarchicalFacialSegmentation(SimilarityMatrix,n_levels,type,runs)
% This function creates sub-modules of the face through iterative runs of
% the spectral clusteinrg (Ng et al.)
% INPUT
% corr_matrix: n x n, m participants, n coordinates of points
% RefScan: 1x1 meshObj

% check input 
if nargin<3,type='Symmetric Laplacian';end
if nargin<4, runs = 50; end

n_cluster = 1;

n = size(SimilarityMatrix,1);


LabelTotal = ones(n_levels+1,n);

h=waitbar(0,'please wait');
for level=1:n_levels
    
    clustercount = 0;
    
    for cluster=1:n_cluster
        
        a = [datestr(now) ' Level ', num2str(level), ' start', ' cluster ', num2str(cluster), ' start'];
        disp(a)
        
        pos = find(LabelTotal(level,:)==cluster);
        
        if LabelTotal(level:end,pos)>0
            
            if length(pos)>1
                A = SimilarityMatrix(pos,pos);
                
                if length(A)==1
                    label=1;
                else
                    label = mySpectralClustering(A,type,runs);
                end
                label(label==1)=clustercount+1;
                label(label==2)=clustercount+2;
                clustercount = clustercount + 2;
                
                LabelTotal(level+1,pos)=label;
            else
                LabelTotal(level+1:end,pos)=-cluster/(level); % to be adjusted
            end
            
            
        end
        
        a = [datestr(now) ' Level ', num2str(level), ' done', ' cluster ', num2str(cluster), ' done'];
        disp(a)
        
    end
    n_cluster = clustercount;
    waitbar(level/n_levels,h);
end
close(h)
end

function [C_opt,cost]=mySpectralClustering(W,type,runs)
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
    n = size(W,1);
    W(eye(n)==1) = 0;
    switch lower(type)
        case 'ratiocute'
            D = sum(W);
            L = diag(D)-W;
            U = eigs(L,2,'smallestabs');
            U = U(:,2);
            cost = 0;
            C_opt = double(U>0)+1;
            return;
        case 'ncute'
            D = sum(W);
            D = diag(D.^(-0.5));
            I = D*W*D;
            U = eigs(I,2);
            U = U(:,2);
            cost = 0;
            C_opt = double(U>0)+1;
            return;
        case 'symmetric laplacian'
            D = sum(W);
            D = diag(D.^(-0.5));
            I = D*W*D;
            U = eigs(I,3);
            U = D*U(:,2:3);
        case 'weiss' % OPTION TAKEN FOR NATURE GENETICS PAPER
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
        otherwise
    end
    cost = +Inf;
    for i=1:runs
        [C,~,c] = myKmeans(U',k);
        disp(num2str(c));
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



