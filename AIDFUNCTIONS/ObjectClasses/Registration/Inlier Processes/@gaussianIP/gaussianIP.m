classdef gaussianIP < IP
    % This is the abstract interface class for Transformation models
    properties
        PdfP = [];% Inlier Process parameters
    end
    properties (Dependent = true)
        Sigma;
        Mu;
    end
    methods %Constructor
        function obj = gaussianIP(varargin)
          obj = obj@IP(varargin{:});
          if nargin > 0
             Input = find(strcmp(varargin, 'Sigma'));if ~isempty(Input), obj.Sigma = varargin{Input+1}; end
          end
        end
    end  
    methods % Special Setting and Getting
        function sigma = get.Sigma(obj)
            sigma = diag(obj.PdfP);
        end
        function obj = set.Sigma(obj,sigma)
            obj.PdfP = diag(sigma);
        end
        function mu = get.Mu(obj)
            mu = zeros(1,obj.nrSM);
        end
    end
end % classdef

