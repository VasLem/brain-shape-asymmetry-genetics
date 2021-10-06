function clusteringResult = hierarchicalClustering(similarityMat, numberOfLevels, useNormSym, eigVectorsNum, seed, eps, runs)
% A: Input Array Landmarks of size numSamples x 3 x numLandmarks
% numberOfLevels: Number of levels of partitioning, defaults to 5
% useNormSym: Whether to use or not symmetric Laplacian Normalization, defaults to true
% eigVectorsNum: Number of Eigen Vectors to keep for the kmeans++ clustering stage,
%     set to 0 for autodetermination based on PCA parallel analysis, defaults to 0
% eps: The minimal amount to consider higher than 0, defautls to 1e-7
%
% Returns:  A tree-like recursive structure (let it be X) , where each node contains the indices that correspond to it,
%   as well as the 'parts' cell array of two X structures, that make it up, if applicable.

if nargin<6
    eps = 1e-7;
end
if nargin < 7, runs = 50; end
if nargin>=5, rng(seed); end
if nargin<4
    eigVectorsNum = 2;
end
if nargin<2
    numberOfLevels = 5;
end
% From here on, we are based on https://en.wikipedia.org/wiki/Spectral_clustering#Algorithms
% and https://nl.mathworks.com/help/stats/spectralcluster.html


if nargin<3
    % There is reference to the paper describing the Normalized Symmetric Laplacian Matrix
    % in the original publication of Prof. Claes, so we are keeping this as a default
    useNormSym = true;
end


clusteringResult = recursivePartition(similarityMat, numberOfLevels, 1:size(similarityMat,1), 1:size(similarityMat,1), eigVectorsNum, useNormSym, eps, runs);
end

function ret = recursivePartition(similarityMat, count, indices, rootIndices, eigVectorsNum, useNormSym, eps, runs)
if isnan(indices)
    indices = 1:size(similarityMat,1);
end
ret.indices = indices;
ret.rootIndices = rootIndices;
ret.parts = {};
if length(indices)<eigVectorsNum
    return
end
if count > 0
    disp(['Level ' num2str(count)]);
    similarityMat(eye(size(similarityMat,1))==1) = 0;
    %1. Compute Laplacian Matrix
    degreeMat = sparse(1:size(similarityMat,1), 1:size(similarityMat,1), sum(similarityMat,1));
    laplacianMat = degreeMat - similarityMat;
    
    %2.Compute Normalized Laplacian Matrix
    if ~useNormSym
        invDegreeMat = 1./max(eps, degreeMat);
        normalizedLaplacianMat = laplacianMat * invDegreeMat ;
        
    else % Compute Normalized Symmetric Laplacian Matrix
        sqrtInvDegreeMat = spdiags(1./max(eps, sqrt(diag(degreeMat))), 0, size(similarityMat,1),  size(similarityMat,1));
        normalizedLaplacianMat = sqrtInvDegreeMat * laplacianMat * sqrtInvDegreeMat;
    end
    %3. Find the eigenvectors of the Normalized Symmetric Laplacian Matrix
    
    sv = randn(size(normalizedLaplacianMat,1),1); % pass random vector to eigs so that to get reproducible results (if using rng from Matlab)
    
    [V,~] = eigs(double(normalizedLaplacianMat),eigVectorsNum, eps, 'StartVector',sv);
    
    if useNormSym
        %Normalize V so that rows sum to 1
        V = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2))); 
    end
    x = V';
    cost = +Inf;
    k = 2;
    for i=1:runs
        [split,~,c] = myKmeans(x,k);
        if c < cost
            s_opt = split;
            cost = c;
        end
    end
    split = s_opt;
    lab1 = find(split == 1);
    lab2 = find(split == 2);
    ret.parts{1} = recursivePartition(similarityMat(lab1,lab1), count-1, lab1, rootIndices(lab1), eigVectorsNum, useNormSym, eps, runs);
    ret.parts{2} = recursivePartition(similarityMat(lab2,lab2), count-1, lab2, rootIndices(lab2), eigVectorsNum, useNormSym, eps, runs);
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

