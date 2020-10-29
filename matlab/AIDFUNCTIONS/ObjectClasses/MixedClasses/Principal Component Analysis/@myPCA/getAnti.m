function out = getAnti(obj,coeff,val)
         if nargin < 3
            val = -1;
         end
         out = multiplyCoeff(obj,coeff,val);
end