function out = getDirection(obj,startcoeff,endcoeff)
% out = getDirection(obj,coeff1,coeff2)
% Determine the Direction between two entities in PCA space
% INPUT
% obj = PCA space object
% startcoeff = PCA coefficients of the first entity or startcoeff = [] being the average entity
% endcoeff = PCA coefficients of the second entity 
% OUTPUT
% out = direction
%
% created by Peter Claes
         if isempty(startcoeff), startcoeff = obj.AvgCoeff; end
         if size(startcoeff,1)==1, startcoeff = startcoeff';end
         if size(endcoeff,1)==1, endcoeff = endcoeff';end
         out = endcoeff-startcoeff;
         out = out/norm(out);
end