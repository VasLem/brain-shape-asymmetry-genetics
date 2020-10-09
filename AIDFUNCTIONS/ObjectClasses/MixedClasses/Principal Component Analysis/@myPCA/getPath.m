function out = getPath(obj,index)
% out = getPath(obj,index)
% Determine the Path to change the vector elements located at index in a PCA space
% A path is a linear combination of eigenvectors that change the value of a
% particular element by 1;
% INPUT
% obj = PCA space object
% index = place of elements in vector representation that need to change in
% if length(index)==1, a single path is returned
% if length(index)>1, multiple perpendicular (independent) paths are returned
% OUTPUT
% out = path
%
% created by Peter Claes
         U = obj.EigVec(index,:);
         A = diag(obj.EigVal);
         out = A*U'/(U*A*U');
end