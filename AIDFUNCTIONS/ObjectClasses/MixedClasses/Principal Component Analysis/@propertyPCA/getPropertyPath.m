function out = getPropertyPath(obj,prop,indprop)
         if isempty(prop)
            prop = (1:obj.nrP);
         elseif ischar(prop)
             prop = strmatch(prop,obj.PropNames);
         else
             %prop = prop;
         end
         if nargin < 3
            indprop = [];
         elseif ischar(indprop)
             indprop = strmatch(indprop,obj.PropNames);
         else
             %indprop = indprop;
         end
         index = [prop setdiff(indprop,prop)]+obj.nrMC;
         out = getPath(obj,index);
         out = out(:,1:length(prop));  
end