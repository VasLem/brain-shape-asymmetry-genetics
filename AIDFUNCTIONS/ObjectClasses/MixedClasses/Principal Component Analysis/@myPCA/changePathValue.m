function out = changePathValue(obj,coeff,delta,path) %#ok<INUSL>
% out = changePathValue(obj,coeff,delta,path)
% created by Peter Claes
    if size(coeff,1) == 1, coeff = coeff'; end
    out = coeff+path*delta;
end