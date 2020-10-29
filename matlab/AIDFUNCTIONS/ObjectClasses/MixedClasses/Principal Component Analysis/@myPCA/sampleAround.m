function out = sampleAround(obj,coeff,dir,D,type)
% out = sampleAround(obj,coeff,dir,D)
% out = sampleAround(obj,coeff,dir,D,type)
% sample around a given entity in PCA space
% INPUT
% obj = PCA space object
% coeff = set of coefficients you want to sample around, coeff = [] sample aroung average
% dir = direction of sampling; dir = [] (random direction); dir = scalar (PC direction) sign(scalar) gives which way; dir = vector (given direction);
% D = distance along the direction to sample
% Type = type of distance( Euclidean or Mahalanobis), default Mahalanobis
% OUTPUT
% out = new set of coeff, resulting from the sampling
%
% Created by Peter Claes
         if isempty(coeff), coeff = obj.AvgCoeff; end
         if size(coeff,1)==1, coeff = coeff'; end
         if nargin < 5, type = 'Mahalanobis';end
         if isempty(dir)% random direction
            dir = -1+2*rand(1,obj.nrEV);dir = (dir/norm(dir))';
         elseif isscalar(dir)% according to a principal component
            pc = dir;
            dir = zeros(1,obj.nrEV)';
            switch sign(pc)
                case 1
                    dir(abs(pc)) = 1;
                case -1
                    dir(abs(pc)) = -1;
            end
         else
             if size(dir,1)==1, dir = dir'; end
         end
         switch lower(type)
             case 'euclidean'
                 A = sum((dir).^2);
             case 'mahalanobis'
                 A = sum((dir./obj.EigStd).^2);
         end
         X = sqrt((D^2)/A);
         out = coeff+dir*X;
end