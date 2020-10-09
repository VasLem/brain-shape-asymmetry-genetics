function out = showShapeEffect(obj,coeff,type,name)
% out = showShapeEffect(obj,coeff,type)
% created by Peter Claes
    vec = Coeff2Vec(obj,coeff);
    if nargin <4, name = 'Shape Effect'; end
    if nargin < 3, type = 'Distance'; end
    out = showShapeEffect(obj.Model,vec(1:obj.nrMC),type,name);
end