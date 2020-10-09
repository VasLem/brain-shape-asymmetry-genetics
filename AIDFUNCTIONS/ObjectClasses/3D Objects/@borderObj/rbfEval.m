function value = rbfEval(obj,input,varargin)
         if isempty(obj.RBF)
            rbfFit(obj,varargin{:});
         end
         value = eval(obj.RBF,input,varargin{:});
         %pl = fastrbf_pointeval(obj.RBF,pl,'messages',0);
         %eval = pl.Value;
end