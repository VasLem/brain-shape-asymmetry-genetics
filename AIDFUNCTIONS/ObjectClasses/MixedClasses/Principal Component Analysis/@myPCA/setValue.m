function out = setValue(obj,coeff,index,value,path)
% out = setValue(obj,coeff,prop,value,path)
% created by Peter Claes
    if isempty(coeff), coeff = obj.AvgCoeff; end
    if isempty(path),path = getPath(obj,index);end
    vec = Coeff2Vec(obj,coeff);
    delta = value-vec(index);
    out = changePathValue(obj,coeff,delta,path);
end