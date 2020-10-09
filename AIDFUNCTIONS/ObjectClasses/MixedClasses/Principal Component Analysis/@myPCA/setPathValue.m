function out = setPathValue(obj,coeff,index,value,path)
% out = setPathValue(obj,coeff,index,value)
% out = setPathValue(obj,coeff,index,value,path)
% creates a new set of coefficients of an entity for which the element
% values are equal to the given value, in other words the entity's element
% value have been change acording to the associated path.
% INPUT
% obj = PCA space object
% ent = starting entity or coeff = [] the average entity
% index = index of the vector elements to change
% value = the desired value of the vector elements
% path = the path to change along, if not given will be computed first
% OUTPUT
% out = new set of coefficients for an entity having the approriate values
% in their vector representation. (see vec = Coeff2Vec(obj,out));
%
% created by Peter Claes
    if isempty(coeff), coeff = obj.AvgVec; end
    if nargin < 5,path = getPath(obj,index);end
    vec = Coeff2Vec(obj,coeff);
    delta = vec(index)-value;
    out = changePathValue(obj,coeff,delta,path);
end