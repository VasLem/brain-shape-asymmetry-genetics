function W = computeSmoothWeightsv2(obj,type,varargin)
% computeWeight - compute a weight matrix
%
%   W = computeWeights(obj,weights,type,normalise);
%
%   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
%   connected in the mesh.
%
%   type is either 
%       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
%       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
%           i and j.
%       'texturedistance': W(i,j) = 1/d_ij^2 where d_ij is distance between RGBvertex
%           i and j.
%       'conformal': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j)
%   
    switch type
        case 'combinatorial'
            n = size(obj.Vertices,2);
            W = obj.Adjacency+speye(n);
            D = spdiags(full(sum(W,2).^(-1)),0,n,n);
            W = D*W;            
        case 'distance'
            W = vertexAdjacency(obj);
            W(W>0) = 1./W(W>0);
            W = (W+W')/2;
            W = diag(sum(W,2).^(-1)) * W;
        case 'texturedistance'
            W = textureAdjacency(obj);
            W(W>0) = 1./W(W>0);
            W = (W+W')/2;
            W = diag(sum(W,2).^(-1)) * W;
        case 'indexeddistance'
            W = indexedAdjacency(obj);
            W(W>0) = 1./W(W>0);
            W = (W+W')/2;
            W = diag(sum(W,2).^(-1)) * W;
        case 'conformal'
            W = angleAdjacency(obj);
            W = diag(sum(W,2).^(-1)) * W;
        case 'functiondistance'
            W = functionAdjacencyv2(obj,varargin{1}');
            W = diag(sum(W,2).^(-1)) * W;
        otherwise
            error('Unknown Type');
    end
end