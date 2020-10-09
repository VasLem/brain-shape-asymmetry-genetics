classdef uniformOP < OP
    % This is the abstract interface class for Transformation models
    properties 
        PdfP = 0;% zero uniform distribution 
    end
    properties (Dependent = true)
        Level;% Level of the unoform distribution
    end
    methods %Constructor
        function obj = uniformOP(varargin)
          obj = obj@OP(varargin{:});
          if nargin>0
              Input = find(strcmp(varargin, 'Level'));if ~isempty(Input), obj.Level = varargin{Input+1}; end
          end
        end
    end  
    methods % Special Setting and Getting
        function level = get.Level(obj)
            level = obj.PdfP;
        end
        function obj = set.Level(obj,level)
            obj.PdfP = level;
        end
    end
end % classdef

