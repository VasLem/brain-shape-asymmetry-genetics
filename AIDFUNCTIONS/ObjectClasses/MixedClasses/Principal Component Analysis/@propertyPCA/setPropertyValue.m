function out = setPropertyValue(obj,coeff,prop,value,path)
% out = setPropertyValue(obj,coeff,prop,value,path)
% created by Peter Claes
    if isempty(coeff), coeff = obj.AvgCoeff; end
    if isempty(path),path = getPropertyPath(obj,prop);end
    vec = Coeff2Vec(obj,coeff);
    delta = value-vec(prop+obj.nrMC);
    out = changePropertyValue(obj,coeff,delta,path);
end