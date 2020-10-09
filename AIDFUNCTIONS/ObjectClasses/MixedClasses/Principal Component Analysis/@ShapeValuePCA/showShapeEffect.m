function out = showShapeEffect(obj,coeff,type,name)
% out = showShapeEffect(obj,coeff,type)
% created by Peter Claes
     if isscalar(coeff)
         pc = coeff;
         coeff = obj.AvgCoeff;
         coeff(pc) = 1;
     end
     vec = Coeff2Vec(obj,coeff);
     if nargin < 3, type = 'Distance'; end
     out = showShapeEffect(obj.Shape,vec(1:obj.nrSC),type,name);
end