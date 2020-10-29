function [out,outPP,outPar,angle] = getDistance(obj,startcoeff,endcoeff,type,Direction)
% out = getDistance(obj,startcoeff,endcoeff,type)
% Determine the Distance between two entities in PCA space
% INPUT
% obj = PCA space object
% startcoeff = PCA coefficients of the first entity or startcoeff = [] being the average entity
% endcoeff = PCA coefficients of the second entity 
% type = 'Mahalanobis' (default) or 'Euclidean'
% OUTPUT
% out = distance
%
% created by Peter Claes
         if nargin < 5, Direction = []; end;
         if nargin < 4, type = 'Mahalanobis'; end
         if isempty(startcoeff),startcoeff = obj.AvgCoeff;end
         if size(startcoeff,1)==1, startcoeff = startcoeff';end      
         if size(endcoeff,1)==1, endcoeff = endcoeff';end
         switch lower(type)
             case 'euclidean'
                out = sqrt(sum((endcoeff-startcoeff).^2));
             case 'mahalanobis'
                out = sqrt(sum(((endcoeff-startcoeff)./obj.EigStd).^2));
         end
         if isempty(Direction), outPP = []; outPar = [];return; end
         CoeffD = (endcoeff-startcoeff);
         CoeffD = CoeffD/norm(CoeffD);
         [angle,~,~] = getAngle(obj,CoeffD,Direction,type);
         outPP = abs(sin(acos(angle))*out);
         outPar = angle*out;
end