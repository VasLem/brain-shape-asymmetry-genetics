classdef dummyOP < OP
    % This is a dummy Outlier class, as if no outliers are defined
    properties
        PdfP=[];% Outlier Process parameters
    end
    methods %Constructor
        function obj = dummyOP(varargin)
          obj = obj@OP(varargin{:});
        end
    end        
    methods %Interfaces Functions
        function out = eval(obj,Tmodel) %#ok<INUSD>
                 if nargout == 1, out = 0; return; end
                 obj.Evaluation = 0;
        end
        function out = derivate(obj,Tmodel) %#ok<INUSD>
                 grad = 0;
                 if nargout == 1, out = grad; return; end
                 obj.Derivative = grad;        
        end
        function out = update(obj,Tmodel) %#ok<INUSD>
                 if nargout == 1, out = []; end
        end
        function out = pdf(obj,Tmodel)
                 val = zeros(1,Tmodel.nrV);
                 if nargout == 1, out = val; return; end
                 obj.PdfEvaluation = val;
        end
    end   
end % classdef

