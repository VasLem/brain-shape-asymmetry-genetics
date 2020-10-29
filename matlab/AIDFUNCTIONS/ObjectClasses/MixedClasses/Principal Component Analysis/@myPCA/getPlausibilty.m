function [out,percout] = getPlausibilty(obj,coeff)
% out = getPlausibilty(obj,coeff)
% Determine the plausibility of an entity in PCA space
% INPUT
% obj = PCA space object
% coeff = PCA coefficients of the entity
% OUTPUT
% out = plausibility
% percout = plausibility expressed in percentage according to maximum
% plausibility
% created by Peter Claes
     if size(coeff,1)==1, coeff = coeff';end
     out = (1/sqrt(2*pi))*exp(-0.5*(sum((coeff./obj.EigStd).^2)));
     maxp = (1/sqrt(2*pi))*exp(-0.5*(sum((zeros(size(coeff))./obj.EigStd).^2)));
     percout = (out/maxp)*100;
end