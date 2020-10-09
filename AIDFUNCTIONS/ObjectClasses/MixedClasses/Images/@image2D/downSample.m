function out = downSample(obj,factor)
         if nargout == 1
            obj = clone(obj);
            out = obj;
         end
         obj.Image = obj.Image(1:factor:end,1:factor:end,:);
end