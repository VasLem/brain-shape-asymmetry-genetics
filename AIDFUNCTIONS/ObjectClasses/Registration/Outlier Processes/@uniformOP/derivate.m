function out = derivate(obj,Tmodel) %#ok<INUSD>
         grad = 0;
         if nargout == 1, out = grad; return; end
         obj.Derivative = grad;
end