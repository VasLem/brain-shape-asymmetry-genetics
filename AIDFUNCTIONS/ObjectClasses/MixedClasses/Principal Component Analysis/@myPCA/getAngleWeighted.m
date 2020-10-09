function [out,out2,percout] = getAngleWeighted(obj,coeff1,coeff2,type,Cor)
% out = getAngle(obj,coeff1,coeff2,type)
% Determine the Angle between two entities in PCA space
% INPUT
% obj = PCA space object
% coeff1 = PCA coefficients of the first entity
% coeff2 = PCA coefficients of the second entity
% type = 'Mahalanobis' (default) or 'Euclidean'
% OUTPUT
% out = angle
% percout = angle expressed in percentage similar
%
% created by Peter Claes and modified by Dorothy Gors
         if nargin < 4, type = 'Mahalanobis'; end 
         if size(coeff1,1)==1, coeff1 = coeff1';end
         if size(coeff2,1)==1, coeff2 = coeff2';end

         switch lower(type)
             case 'euclidean'
                % do nothing
             case 'mahalanobis'
                coeff1 = coeff1./obj.EigStd;
                coeff2 = coeff2./obj.EigStd;
             case 'weighted'
                coeff1 = (coeff1./obj.EigStd).*abs(Cor');
                coeff2 = (coeff2./obj.EigStd).*abs(Cor');
         end
         T = coeff1'*coeff2;
         N = sqrt((coeff1'*coeff1)*(coeff2'*coeff2));
         out = T/N;
         out2 = acosd(out);
         percout = ((out+1)/2)*100;
end