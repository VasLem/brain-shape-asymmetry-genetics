function coeff = sampleRandom(obj,type,v)
% coeff = sampleRandom(obj,type,v)
% randomly sampling space coeff
% INPUT
% obj = PCA object
% type = 'Uniform' (Default) or 'Normal'
% v = factor to push (< 1) closer or further (> 1) away from average
% OUTPUT
% coeff = random set of coeff;
% Created by Peter Claes
         if nargin < 2, type = 'Uniform'; v = 1;end
         if nargin < 3, v = 1;end
         switch type
             case 'Normal'
                 coeff = randn(1,obj.nrEV);
             case 'Uniform'
                 coeff = -1+2*rand(1,obj.nrEV);
         end      
         coeff = v*coeff.*obj.EigStd';
end