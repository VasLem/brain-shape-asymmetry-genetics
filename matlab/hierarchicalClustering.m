function clusteringResult = hierarchicalClustering(A, numberOfLevels, useNormSym, eigVectorsNum, seed, eps)
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
if nargin>=5
    rng(seed);
end
if nargin<4
    eigVectorsNum = 0;
end
if nargin<2
    numberOfLevels = 5;
end
% Steps for hierarchical Clustering
% Quasi Landmarks have to be already downsampled and aligned to each other
% 1.Compute the pairwise 3D correlation matrix
meanConditionalOnLandmarks = mean(A, 3);
stdConditionalOnLandmarks = std(A, 0, 3);

B = (A - meanConditionalOnLandmarks)./max(eps, stdConditionalOnLandmarks);
resB = reshape(B,[size(A,1), size(A,2) * size(A,3)]);
similarityMat = abs(resB*resB');
% From here on, we are based on https://en.wikipedia.org/wiki/Spectral_clustering#Algorithms 
% and https://nl.mathworks.com/help/stats/spectralcluster.html

%2. Compute Laplacian Matrix
if nargin<3
    % There is reference to the paper describing the Normalized Symmetric Laplacian Matrix 
    % in the original publication of Prof. Claes, so we are keeping this as a default
    useNormSym = true;
end
degreeMat = diag(sum(similarityMat,1));
laplacianMat = degreeMat - similarityMat;

%3.Compute Normalized Laplacian Matrix
if ~useNormSym 
    invDegreeMat = diag(1./max(eps, diag(degreeMat)));
    normalizedLaplacianMat = laplacianMat * invDegreeMat ;
    
else % Compute Normalized Symmetric Laplacian Matrix
    sqrtInvDegreeMat = diag(1./max(eps, sqrt(diag(degreeMat))));
    normalizedLaplacianMat = sqrtInvDegreeMat * laplacianMat * sqrtInvDegreeMat;
end
%4. Find the eigenvectors of the Normalized Symmetric Laplacian Matrix


if eigVectorsNum == 0
    %4a Find optimal number of components using parallelAnalysis
    [latent, ~, latentHigh] = parallelAnalysis(normalizedLaplacianMat);
    eigVectorsNum = sum(latent > latentHigh');
end
sv = randn(size(normalizedLaplacianMat,1),1); % pass random vector to eigs so that to get reproducible results (if using rng from Matlab)

[V,~] = eigs(double(normalizedLaplacianMat),eigVectorsNum,1e-10,'StartVector',sv);

if useNormSym
    %Normalize V so that rows sum to 1
    vnorm = vecnorm(V,2,2); % Use vecnorm for improved numerical precision
    nonZeroVNorm = vnorm > 0;
    vnorm(~nonZeroVNorm) = 0;
    if any(nonZeroVNorm)
        V(nonZeroVNorm,:) = V(nonZeroVNorm,:)./vnorm(nonZeroVNorm);
    end
end
clusteringResult = recursivePartition(V', numberOfLevels, 1:size(V,1));
end

function ret = recursivePartition(x, count, indices)
    if isnan(indices)
        indices = 1:size(x,2);
    end
    ret.indices = indices;
    ret.parts = {};
    if count > 0
        [split,~] = kmeans(x,2);
        if any(split == 2)
            if sum(split == 1)>sum(split == 2)
                c1 = 1;
                c2 = 2;
            else
                c1 = 2;
                c2 = 1;
            end
            ret.parts{c1}= recursivePartition(x(split == 1), count-1, indices(split==1));
            ret.parts{c2} = recursivePartition(x(split==2), count-1, indices(split==2));
        end
    end
end
   


function [L,C] = kmeans(X,k)
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
% _CITE: Laurent S (2021). k-means++ (https://www.mathworks.com/matlabcentral/fileexchange/28804-k-means), MATLAB Central File Exchange. Retrieved September 8, 2021.
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
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
    
end
end

