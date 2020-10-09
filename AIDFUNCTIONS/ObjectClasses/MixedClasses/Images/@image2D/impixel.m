function out = impixel(obj,Field)
         if nargin == 1
             im = obj.Image;
         else
             im = obj.(Field);
         end
         out = impixel(im); 
end